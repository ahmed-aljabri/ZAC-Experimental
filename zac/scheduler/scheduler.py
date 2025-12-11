import time
import math
import rustworkx as rx

class Scheduler_mixin:
    def scheduling(self):
        """
        solve gate scheduling problem for all-commutable gate cases by graph coloring algorithm
        """
        t_s = time.time()
        self.gate_scheduling = []
        if self.has_dependency:
            if self.scheduling_strategy == "asap":
                result_scheduling = self.asap()
            else:
                result_scheduling = [[i] for i in range(len(self.g_q))]
        else:
            result_scheduling = self.graph_coloring()
        
        # capacity and pair capacity
        max_gate_num = 0
        for zone in self.architecture.entanglement_zone:
            slm = self.architecture.dict_SLM[zone[0]]
            max_gate_num += (slm.n_r * slm.n_c)
        pair_capacity = max_gate_num // 2 if max_gate_num > 0 else 0
        
        self.gate_scheduling_idx = []
        estimated_qubit_loc = None
        last_ryd_site = None
        zone_label = None
        if getattr(self, "locality_reorder", False) or getattr(self, "locality_cluster_split", False):
            estimated_qubit_loc = self._estimate_qubit_locations_for_reorder()
            last_ryd_site = self._init_last_ryd_site(estimated_qubit_loc)
            zone_label = self._zone_labels()

        for stage_idx, gates in enumerate(result_scheduling):
            layer_gates = gates
            if getattr(self, "locality_reorder", False) or getattr(self, "locality_cluster_split", False):
                layer_gates = self._reorder_layer_by_locality(stage_idx, layer_gates, estimated_qubit_loc, last_ryd_site, zone_label)
                self._update_estimates_after_layer(layer_gates, estimated_qubit_loc, last_ryd_site)
            if len(layer_gates) < max_gate_num or (not getattr(self, "locality_cluster_split", False)) or pair_capacity == 0:
                self.gate_scheduling_idx.append(layer_gates)
            else:
                num_layer = math.ceil(len(layer_gates) / pair_capacity)
                cluster_size = math.ceil(len(layer_gates) / num_layer)
                clusters = [layer_gates[i:i+cluster_size] for i in range(0, len(layer_gates), cluster_size)]
                for cluster in clusters:
                    self.gate_scheduling_idx.append(cluster)
                self._trace_layer(stage_idx, layer_gates, zone_label, pair_capacity, clusters, last_ryd_site)
                continue

        self.gate_scheduling = []
        self.gate_1q_scheduling = []

        for gates in self.gate_scheduling_idx:
            tmp = [self.g_q[i] for i in gates]
            self.gate_scheduling.append(tmp)
            self.gate_1q_scheduling.append([])
            for gate_idx in gates:
                if gate_idx in self.dict_g_1q_parent:
                    for gate_1q in self.dict_g_1q_parent[gate_idx]:
                        self.gate_1q_scheduling[-1].append(gate_1q)
        
        self.runtime_analysis["scheduling"] = time.time()- t_s
        exceeding = [len(gates) for gates in self.gate_scheduling_idx if len(gates) > max_gate_num]
        print(f"[TRACE] scheduling: rydberg_layers={len(self.gate_scheduling_idx)}, capacity_per_layer={max_gate_num}, layers_exceeding={len(exceeding)}")
        if exceeding:
            print(f"[TRACE] scheduling: lengths_exceeding_cap={exceeding[:10]}")
        print("[INFO]               Time for scheduling: {}s".format(self.runtime_analysis["scheduling"]))

    def asap(self):
        # as soon as possible algorithm for self.g_q
        gate_scheduling = []
        list_qubit_time = [0 for i in range(self.n_q)]
        for i, gate in enumerate(self.g_q):
            tq0 = list_qubit_time[gate[0]]
            tq1 = list_qubit_time[gate[1]]
            tg = max(tq0, tq1)
            if tg >= len(gate_scheduling):
                gate_scheduling.append([])
            gate_scheduling[tg].append(i)

            tg += 1
            list_qubit_time[gate[0]] = tg
            list_qubit_time[gate[1]] = tg
        return gate_scheduling
    
    def graph_coloring(self):
        graph = rx.PyGraph()
        graph.add_nodes_from(list(range(self.n_q)))
        for edge in self.g_q:
            graph.add_edge(edge[0], edge[1], edge)
        edge_colors = rx.graph_misra_gries_edge_color(graph)
        max_color = 0
        for i in edge_colors:
            max_color = max(max_color, edge_colors[i])
        gate_scheduling = [[] for i in range(max_color + 1)]
        for i in range(len(self.g_q)):
            gate_scheduling[edge_colors[i]].append(i)
        return gate_scheduling

    def _estimate_qubit_locations_for_reorder(self):
        locs = [None for _ in range(self.n_q)]
        slm_ids = list(self.architecture.storage_zone)
        idx_slm = 0
        r = 0
        c = 0
        for q in range(self.n_q):
            slm = self.architecture.dict_SLM[slm_ids[idx_slm]]
            locs[q] = (slm.idx, r, c)
            c += 1
            if c >= slm.n_c:
                c = 0
                r += 1
            if r >= slm.n_r and idx_slm + 1 < len(slm_ids):
                idx_slm += 1
                slm = self.architecture.dict_SLM[slm_ids[idx_slm]]
                r = 0
                c = 0
        for q in range(self.n_q):
            if locs[q] is None:
                locs[q] = locs[q % len(locs)]
        return locs

    def _init_last_ryd_site(self, locs):
        last = [None for _ in range(self.n_q)]
        for q, loc in enumerate(locs):
            cands = self.architecture.nearest_entanglement_site(loc[0], loc[1], loc[2], loc[0], loc[1], loc[2])
            if len(cands) > 0:
                last[q] = cands[0]
        return last

    def _reorder_layer_by_locality(self, stage_idx, layer_gates, locs, last_sites, zone_label):
        def storage_loc(loc):
            slm = self.architecture.dict_SLM[loc[0]]
            if slm.entanglement_id == -1:
                return loc
            return self.architecture.nearest_storage_site(loc[0], loc[1], loc[2])

        def gate_key(idx_gate):
            q0, q1 = self.g_q[idx_gate]
            site0 = last_sites[q0]
            site1 = last_sites[q1]
            if site0 is None:
                loc0 = storage_loc(locs[q0])
                c0 = self.architecture.nearest_entanglement_site(loc0[0], loc0[1], loc0[2], loc0[0], loc0[1], loc0[2])
                site0 = c0[0] if len(c0) > 0 else None
            if site1 is None:
                loc1 = storage_loc(locs[q1])
                c1 = self.architecture.nearest_entanglement_site(loc1[0], loc1[1], loc1[2], loc1[0], loc1[1], loc1[2])
                site1 = c1[0] if len(c1) > 0 else None
            zone = "unknown"
            if site0 is not None:
                zid = self.architecture.dict_SLM[site0[0]].entanglement_id
                zone = zone_label.get(zid, "unknown")
            elif site1 is not None:
                zid = self.architecture.dict_SLM[site1[0]].entanglement_id
                zone = zone_label.get(zid, "unknown")
            span = 0
            if site0 is not None and site1 is not None:
                span = self.architecture.distance(site0[0], site0[1], site0[2], site1[0], site1[1], site1[2])
            row_proxy = min([site0[1] if site0 else 1e9, site1[1] if site1 else 1e9])
            zone_order = {"top": 0, "bot": 1, "unknown": 2} if stage_idx % 2 == 0 else {"bot": 0, "top": 1, "unknown": 2}
            zone_pri = zone_order.get(zone, 2)
            return (zone_pri, row_proxy, span, idx_gate)

        if all(last_sites[self.g_q[idx][0]] is None and last_sites[self.g_q[idx][1]] is None for idx in layer_gates):
            return layer_gates
        return sorted(layer_gates, key=gate_key)

    def _update_estimates_after_layer(self, layer_gates, locs, last_sites):
        def storage_loc(loc):
            slm = self.architecture.dict_SLM[loc[0]]
            if slm.entanglement_id == -1:
                return loc
            return self.architecture.nearest_storage_site(loc[0], loc[1], loc[2])
        for idx_gate in layer_gates:
            q0, q1 = self.g_q[idx_gate]
            loc0 = storage_loc(locs[q0])
            loc1 = storage_loc(locs[q1])
            cands = self.architecture.nearest_entanglement_site(loc0[0], loc0[1], loc0[2], loc1[0], loc1[1], loc1[2])
            if len(cands) == 0:
                continue
            site = cands[0]
            last_sites[q0] = site
            last_sites[q1] = (site[0]+1, site[1], site[2])
            locs[q0] = last_sites[q0]
            locs[q1] = last_sites[q1]

    def _zone_labels(self):
        y_zone = []
        for zone in self.architecture.entanglement_zone:
            slm = self.architecture.dict_SLM[zone[0]]
            y_zone.append((slm.location[1], self.architecture.dict_SLM[zone[0]].entanglement_id))
        y_zone = sorted(y_zone, key=lambda x: x[0])
        label = dict()
        if len(y_zone) > 0:
            label[y_zone[0][1]] = "top"
        if len(y_zone) > 1:
            label[y_zone[1][1]] = "bot"
        for _, zid in y_zone[2:]:
            label[zid] = "unknown"
        return label

    def _trace_layer(self, stage_idx, layer_gates, zone_label, pair_capacity, clusters, last_sites):
        if zone_label is None or last_sites is None:
            return
        keys = []
        zone_counts = {"top": 0, "bot": 0, "unknown": 0}
        row_vals = []
        for idx_gate in layer_gates:
            q0, q1 = self.g_q[idx_gate]
            site0 = last_sites[q0]
            site1 = last_sites[q1]
            zone = "unknown"
            if site0 is not None:
                zid = self.architecture.dict_SLM[site0[0]].entanglement_id
                zone = zone_label.get(zid, "unknown")
            elif site1 is not None:
                zid = self.architecture.dict_SLM[site1[0]].entanglement_id
                zone = zone_label.get(zid, "unknown")
            zone_counts[zone] = zone_counts.get(zone, 0) + 1
            row_proxy = min([site0[1] if site0 else 1e9, site1[1] if site1 else 1e9])
            row_vals.append(row_proxy)
            span = 0
            if site0 is not None and site1 is not None:
                span = self.architecture.distance(site0[0], site0[1], site0[2], site1[0], site1[1], site1[2])
            keys.append((zone, row_proxy, span, idx_gate))
        print(f"[TRACE] stage {stage_idx}: gates={len(layer_gates)}, pair_cap={pair_capacity}, clusters={[len(c) for c in clusters]}")
        total = len(layer_gates) if len(layer_gates) > 0 else 1
        print(f"[TRACE]   zone mix: top={zone_counts['top']/total:.2f}, bot={zone_counts['bot']/total:.2f}, unknown={zone_counts['unknown']/total:.2f}, avg_row={sum(row_vals)/total if row_vals else 0:.2f}")
        if keys:
            print(f"[TRACE]   first_keys={keys[:10]}")
