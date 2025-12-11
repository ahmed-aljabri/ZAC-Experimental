import time
import subprocess
import heapq
import math
from random import seed, randrange
from copy import deepcopy
from itertools import chain
# from memory_profiler import profile

seed(0)

class Router_mixin:

    # constant, the distance of AOD row and col to some trap. We use 1um here.
    PARKING_DIST = 1
    
    # @profile
    def route_qubit(self):
        """
        generate rearrangement layers between two Rydberg layers
        """
        self._init_aod_heaps()
        self.rydberg_dependency = [0 for i in range(len(self.architecture.entanglement_zone))]
        time_mis = 0
        self.qubit_dependency = [0 for i in range(self.n_q)]
        self.site_dependency = dict()
        self.write_initial_instruction()
        
        

        for layer in range(len(self.gate_scheduling)):
            # extract sets of movement that can be perform simultaneously
            t_s = time.time()
            self.route_qubit_mis(layer)
            time_mis += (time.time() - t_s)
            print("[INFO] ZAC: Solve for Rydberg stage {}/{}. mis time={:2f}".format(layer+1, len(self.gate_scheduling), time_mis))
        self.flatten_rearrangment_instruction()
        self.runtime_analysis["routing"] = time_mis

    def route_qubit_mis(self, layer:int):
        """
        process layers of movement from storage zone to Rydberg and back to storage zone
        """
        initial_mapping = self.qubit_mapping[2 * layer]
        gate_mapping = self.qubit_mapping[2 * layer+1]
        if layer + 2 < len(self.qubit_mapping):
            final_mapping = self.qubit_mapping[2 * layer+2]
        else:
            final_mapping = None

        remain_graph = [] # consist qubits to be moved
        for gate in self.gate_scheduling[layer]:
            for q in gate:
                if initial_mapping[q] != gate_mapping[q]:
                    remain_graph.append(q)
        
        if not(self.routing_strategy == "mis" or self.routing_strategy == "maximalis"):
            remain_graph = sorted(remain_graph, key=lambda x: (math.dist(self.architecture.exact_SLM_location_tuple(initial_mapping[x]),\
                                                                         self.architecture.exact_SLM_location_tuple(gate_mapping[x]))), reverse=True)
        
        id_layer_start = len(self.result_json['instructions'])
        batch = 0
        batch = self._process_movement_batches(remain_graph, initial_mapping, gate_mapping, batch)

        # append a layer for gate execution 
        self.process_gate_layer(layer, gate_mapping)
        # move qubit back to the final location
        if final_mapping is not None:
            if self.dynamic_placement or self.reuse:
                remain_graph = [] # consist qubits to be moved
                for gate in self.gate_scheduling[layer]:
                    for q in gate:
                        if final_mapping[q] != gate_mapping[q]:
                            remain_graph.append(q)
                
                if not(self.routing_strategy == "mis" or self.routing_strategy == "maximalis"):
                    remain_graph = sorted(remain_graph, key=lambda x: (math.dist(self.architecture.exact_SLM_location_tuple(final_mapping[x]), \
                                                                                self.architecture.exact_SLM_location_tuple(gate_mapping[x]))), reverse=True)
                if getattr(self, "interleave_reverse", False):
                    stats = []
                    batch = self._process_movement_batches(remain_graph, gate_mapping, final_mapping, batch, stats)
                    avg_moves = (sum(stats) / len(stats)) if stats else 0
                    print(f"[TRACE] reverse stage {layer}: batches={len(stats)}, avg_moves_per_batch={avg_moves:.2f}")
                else:
                    batch = self._process_movement_batches(remain_graph, gate_mapping, final_mapping, batch)
            else:
                # construct reverse layer
                self.construct_reverse_layer(id_layer_start, gate_mapping, final_mapping)
            self.aod_assignment(id_layer_start)
                
    
    def _process_movement_batches(self, remain_graph: list, initial_mapping: list, final_mapping: list, batch_start: int, stats: list = None):
        """
        Run routing batches either in a single mixed set or grouped by entanglement zone.
        """
        if not remain_graph:
            return batch_start

        if not getattr(self, "zone_batched_routing", False):
            return self._process_single_batch_group(remain_graph, initial_mapping, final_mapping, batch_start, stats)

        zone_to_qubits = dict()
        for q in remain_graph:
            zone_id = self._movement_zone_id(initial_mapping[q], final_mapping[q])
            if zone_id not in zone_to_qubits:
                zone_to_qubits[zone_id] = []
            zone_to_qubits[zone_id].append(q)

        batch = batch_start
        for zone_id in sorted(zone_to_qubits.keys()):
            batch = self._process_single_batch_group(zone_to_qubits[zone_id], initial_mapping, final_mapping, batch, stats)
        return batch

    def _process_single_batch_group(self, remain_graph: list, initial_mapping: list, final_mapping: list, batch_start: int, stats: list = None):
        """
        Core MIS/maximal-IS batching loop for a list of movements.
        """
        batch = batch_start
        pending_qubits = list(remain_graph)
        while pending_qubits:
            vectors = self.graph_construction(pending_qubits, initial_mapping, final_mapping)
            violations = self.collect_violation(vectors)
            if self.routing_strategy == "mis":
                moved_qubits = self.kamis_solve(len(vectors), violations, batch)
            else:
                moved_qubits = self.maximalis_solve(len(vectors), violations)

            set_aod = {pending_qubits[i] for i in moved_qubits} # use to record aods per movement layer
            self.process_movement_layer(set_aod, initial_mapping, final_mapping)
            if stats is not None:
                stats.append(len(set_aod))
            pending_qubits = [q for q in pending_qubits if q not in set_aod]
            batch += 1
        return batch

    def _movement_zone_id(self, start_loc: tuple, target_loc: tuple):
        """
        Choose the entanglement zone driving a move. Prefer the target zone when present,
        otherwise fall back to the source (for reverse moves back to storage).
        """
        zone_to = self.architecture.dict_SLM[target_loc[0]].entanglement_id
        if zone_to != -1:
            return zone_to
        return self.architecture.dict_SLM[start_loc[0]].entanglement_id
    
    def graph_construction(self, remain_graph: list, initial_mapping: list, final_mapping: list):
        vectors = []
        if self.use_window:
            vector_length = min(self.window_size, len(remain_graph))
        else:
            vector_length = len(remain_graph)
        vectors = [(0,0,0,0, ) for _ in range(vector_length)]
        for i, q in enumerate(remain_graph):
            (q_x, q_y) = self.architecture.exact_SLM_location_tuple(initial_mapping[q])
            (site_x, site_y) = self.architecture.exact_SLM_location_tuple(final_mapping[q])
            vectors[i] = (q_x, site_x, q_y, site_y, )
        return vectors
    
    def collect_violation(self, vectors: list):
        violations = []
        for i in range(len(vectors)):
            for j in range(i+1, len(vectors)):
                if not self.compatible_2D(vectors[i], vectors[j]):
                    violations.append((i, j))
        return violations

    def maximalis_solve(self, n, edges):
        """
        solve maximal independent set
        """
        is_node_conflict = [False for _ in range(n)]
        node_neighbors = {i: [] for i in range(n)}
        for edge in edges:
            node_neighbors[edge[0]].append(edge[1])
            node_neighbors[edge[1]].append(edge[0])
        result = []
        for i in range(len(is_node_conflict)):
            if is_node_conflict[i]:
                continue
            else:
                result.append(i)
                for j in node_neighbors[i]:
                    is_node_conflict[j] = True
        return result

    def kamis_solve(self, n, edges, batch):
        """
        solve maximum independent set by KaMIS
        """
        adj = [[] for _ in range(n)]
        for edge in edges:
            adj[edge[0]].append(edge[1])
            adj[edge[1]].append(edge[0])
        for i in range(n):
            adj[i].sort()
        with open(f"mis/tmp/mis_{batch}.in", "w") as f:
            f.write(f"{n} {len(edges)}\n")
            for i in range(n):
                tmp = ""
                for j in adj[i]:
                    tmp += str(j+1)
                    tmp += " "
                tmp += "\n"
                f.write(tmp)
        with open(f"mis/tmp/mis_{batch}.log", "w") as f:
            subprocess.run(
                ["mis/redumis", f"mis/tmp/mis_{batch}.in", "--output", f"mis/tmp/mis_{batch}.out", "--time_limit", "3600"], stdout=f)
        with open(f"mis/tmp/mis_{batch}.out", "r") as f:
            lines = f.readlines()
        return [i for i, line in enumerate(lines) if line.startswith("1")]

    def compatible_2D(self, a, b):
        """
        check if move a and b can be performed simultaneously
        """
        if a[0] == b[0] and a[1] != b[1]:
            return False
        if a[1] == b[1] and a[0] != b[0]:
            return False
        if a[0] < b[0] and a[1] >= b[1]:
            return False
        if a[0] > b[0] and a[1] <= b[1]:
            return False

        if a[2] == b[2] and a[3] != b[3]:
            return False
        if a[3] == b[3] and a[2] != b[2]:
            return False
        if a[2] < b[2] and a[3] >= b[3]:
            return False
        if a[2] > b[2] and a[3] <= b[3]:
            return False
        return True
    

    def write_initial_instruction(self):
        self.result_json['instructions'].clear()
        self.result_json['instructions'].append(
            {
                "type": "init",
                "id": 0,
                "begin_time": 0,
                "end_time": 0,
                "init_locs": [ [i, self.qubit_mapping[0][i][0], self.qubit_mapping[0][i][1], self.qubit_mapping[0][i][2]]
                 for i in range(self.n_q)]
            }
        )

        # process single-qubit gates
        set_qubit_dependency = set()
        inst_idx = len(self.result_json['instructions'])
        list_1q_gate = [gate_1q for gate_1q in self.dict_g_1q_parent[-1]]
        result_gate = []
        for gate_info in list_1q_gate:
            set_qubit_dependency.add(self.qubit_dependency[gate_info[1]])
            self.qubit_dependency[gate_info[1]] = inst_idx
            result_gate.append({
                "name": gate_info[0],
                "q": gate_info[1]
            })
        dependency = { "qubit": []}
        dependency["qubit"] = list(set_qubit_dependency)
        if len(result_gate) > 0:
            self.write_1q_gate_instruction(inst_idx, result_gate, dependency, self.qubit_mapping[0])
            self.result_json['instructions'][-1]["begin_time"] = 0
            self.result_json['instructions'][-1]["end_time"] = (self.architecture.time_1qGate * len(result_gate)) # due to sequential execution

    def process_movement_layer(self, set_aod_qubit: set, initial_mapping: list, final_mapping: list):
        """
        generate layers for row-by-row based atom transfer
        """
        pickup_dict = dict() # key: array and row, value: a list of qubit in the same row
        for q in set_aod_qubit:
            x, y = self.architecture.exact_SLM_location_tuple(initial_mapping[q])
            if y in pickup_dict:
                pickup_dict[y].append(q)
            else:
                pickup_dict[y] = [q]
        list_aod_qubits = []
        list_end_location = []
        list_begin_location = []
        dependency = {
            "qubit": [],
            "site": [],
        }
        inst_idx = len(self.result_json['instructions'])

        set_qubit_dependency = set()
        set_site_dependency = set()
        for dict_key in pickup_dict:
            list_aod_qubits.append(pickup_dict[dict_key])
            row_begin_location = []
            row_end_location = []
            for q in pickup_dict[dict_key]:
                row_begin_location.append([q, initial_mapping[q][0], initial_mapping[q][1], initial_mapping[q][2]])
                row_end_location.append([q, final_mapping[q][0], final_mapping[q][1], final_mapping[q][2]])
                site_key = (final_mapping[q][0], final_mapping[q][1], final_mapping[q][2])
                if site_key in self.site_dependency:
                    set_site_dependency.add(self.site_dependency[site_key])
                site_key = (initial_mapping[q][0], initial_mapping[q][1], initial_mapping[q][2])
                self.site_dependency[site_key] = inst_idx
                set_qubit_dependency.add(self.qubit_dependency[q])
                self.qubit_dependency[q] = inst_idx
            list_begin_location.append(row_begin_location)
            list_end_location.append(row_end_location)
        dependency["qubit"] = list(set_qubit_dependency)
        dependency["site"] = list(set_site_dependency)
        target_zone = None
        for locs in list_end_location:
            for _, a, r, c in locs:
                ent_id = self.architecture.dict_SLM[a].entanglement_id
                if ent_id != -1:
                    target_zone = ent_id
                    break
            if target_zone is not None:
                break
        self.write_rearrangement_instruction(inst_idx, list_aod_qubits, list_begin_location, list_end_location, dependency, target_zone)
    
    def write_rearrangement_instruction(self, inst_idx: int, aod_qubits: list, begin_location: list, end_location: list, dependency: dict, target_zone=None):    
        inst = {
                "type": "rearrangeJob",
                "id": inst_idx,
                "aod_id": -1,
                "target_zone": target_zone,
                "aod_qubits": aod_qubits,
                "begin_locs": begin_location,
                "end_locs": end_location,
                "dependency": dependency
            }
        inst["insts"] = self.expand_arrangement(inst)
        self.result_json['instructions'].append(inst)
    
    def flatten_rearrangment_instruction(self):    
        for inst in self.result_json['instructions']:
            if inst["type"] == "rearrangeJob":
                inst["aod_qubits"] = list(chain.from_iterable(inst["aod_qubits"]))
                inst["begin_locs"] = list(chain.from_iterable(inst["begin_locs"]))
                inst["end_locs"] = list(chain.from_iterable(inst["end_locs"]))
    
    def process_gate_layer(self, layer: int, gate_mapping: list):
        """
        generate a layer for gate execution
        """
        list_gate_idx = self.gate_scheduling_idx[layer]
        list_gate = self.gate_scheduling[layer]
        list_1q_gate = self.gate_1q_scheduling[layer]
        dict_gate_zone = dict()
        for i in range(len(list_gate)):
            slm_idx = gate_mapping[list_gate[i][0]][0]
            zone_idx = self.architecture.dict_SLM[slm_idx].entanglement_id
            if zone_idx not in dict_gate_zone:
                dict_gate_zone[zone_idx] = [i]
            else:
                dict_gate_zone[zone_idx].append(i)
        for rydberg_idx in dict_gate_zone:
            result_gate = [{"id": list_gate_idx[i], "q0": list_gate[i][0], "q1": list_gate[i][1]} for i in dict_gate_zone[rydberg_idx]]
            set_qubit_dependency = set()
            inst_idx = len(self.result_json['instructions'])
            for gate_idx in dict_gate_zone[rydberg_idx]:
                gate = list_gate[gate_idx]
                set_qubit_dependency.add(self.qubit_dependency[gate[0]])
                self.qubit_dependency[gate[0]] = inst_idx
                set_qubit_dependency.add(self.qubit_dependency[gate[1]])
                self.qubit_dependency[gate[1]] = inst_idx
            dependency = { "qubit": [], "rydberg": self.rydberg_dependency[rydberg_idx]}
            self.rydberg_dependency[rydberg_idx] = inst_idx
            dependency["qubit"] = list(set_qubit_dependency)
            self.write_gate_instruction(inst_idx, rydberg_idx, result_gate, dependency)
        
        # process single-qubit gates
        inst_idx = len(self.result_json['instructions'])
        result_gate = []
        set_qubit_dependency = set()
        for gate_info in list_1q_gate:
            set_qubit_dependency.add(self.qubit_dependency[gate_info[1]])
            self.qubit_dependency[gate_info[1]] = inst_idx
            result_gate.append({
                "name": gate_info[0],
                "q": gate_info[1]
            })
        dependency = { "qubit": []}
        dependency["qubit"] = list(set_qubit_dependency)
        if len(result_gate) > 0:
            self.write_1q_gate_instruction(inst_idx, result_gate, dependency, gate_mapping)
    
    def write_gate_instruction(self, inst_idx: int, rydberg_idx: int, result_gate: list, dependency: dict):
        self.result_json['instructions'].append(
            {
                "type": "rydberg",
                "id": inst_idx,
                "zone_id": rydberg_idx,
                "gates": result_gate,
                "dependency": dependency
            }
        )
    
    def write_1q_gate_instruction(self, inst_idx: int, result_gate: list, dependency: dict, gate_mapping: list):
        locs = []
        for gate in result_gate:
            locs.append((gate["q"], gate_mapping[gate["q"]][0], gate_mapping[gate["q"]][1], gate_mapping[gate["q"]][2]))

        self.result_json['instructions'].append(
            {
                "type": "1qGate",
                "unitary": "u3",
                "id": inst_idx,
                "locs": locs,
                "gates": result_gate,
                "dependency": dependency
            }
        )

    def construct_reverse_layer(self, id_layer_start: int, initial_mapping: list, final_mapping: list):
        """
        construct reverse movement layer by processing the forward movement
        """
        id_layer_end = len(self.result_json['instructions'])
        for layer in range(id_layer_start, id_layer_end):
            if self.result_json['instructions'][layer]["type"] == "rydberg":
                break
            else:
                inst_idx = len(self.result_json['instructions'])
                dependency = {
                    "qubit": [],
                    "site": [],
                }
                set_qubit_dependency = set()
                set_site_dependency = set()
                list_aod_qubits = self.result_json['instructions'][layer]["aod_qubits"]
                list_end_location = []
                list_begin_location = []
                for sub_list_qubits in list_aod_qubits:
                    row_begin_location = []
                    row_end_location = []
                    for q in sub_list_qubits:
                        row_begin_location.append([q, initial_mapping[q][0], initial_mapping[q][1], initial_mapping[q][2]])
                        row_end_location.append([q, final_mapping[q][0], final_mapping[q][1], final_mapping[q][2]])
                        site_key = (final_mapping[q][0], final_mapping[q][1], final_mapping[q][2])
                        if site_key in self.site_dependency:
                            set_site_dependency.add(self.site_dependency[site_key])
                        site_key = (initial_mapping[q][0], initial_mapping[q][1], initial_mapping[q][2])
                        self.site_dependency[site_key] = inst_idx
                        set_qubit_dependency.add(self.qubit_dependency[q])
                        self.qubit_dependency[q] = inst_idx

                    list_begin_location.append(row_begin_location)
                    list_end_location.append(row_end_location)
                dependency["qubit"] = list(set_qubit_dependency)
                dependency["site"] = list(set_site_dependency)
                self.write_rearrangement_instruction(inst_idx, 
                                                     list_aod_qubits,
                                                     list_begin_location,
                                                     list_end_location, 
                                                     dependency,
                                                     None)
                

    def aod_assignment(self, id_layer_start: int):
        """
        processs the aod assignment between two Rydberg stages
        """
        list_instruction_duration = [[], []]
        id_layer_end = len(self.result_json['instructions'])
        duration_idx = 0
        list_gate_layer_idx = []
        for idx in range(id_layer_start, id_layer_end):
            if self.result_json['instructions'][idx]["type"] != "rearrangeJob":
                duration_idx = 1
                list_gate_layer_idx.append(idx)
                continue
            duration = self.get_duration(self.result_json['instructions'][idx])
            list_instruction_duration[duration_idx].append((duration, idx))
        list_instruction_duration[0] = sorted(list_instruction_duration[0], reverse=True)
        list_instruction_duration[1] = sorted(list_instruction_duration[1], reverse=True)
        for i in range(2):
            for item in list_instruction_duration[i]:
                duration = item[0]
                inst = self.result_json['instructions'][item[1]]
                begin_time, aod_id = self._pop_aod(inst)
                begin_time = max(begin_time, self.get_begin_time(item[1], inst["dependency"]))
                end_time = begin_time + duration
                inst["dependency"]["aod"] = self.aod_dependency[aod_id]
                self.aod_dependency[aod_id] = item[1]
                inst["begin_time"] = begin_time
                inst["end_time"] = end_time
                inst["aod_id"] = aod_id
                self._push_aod(aod_id, end_time)
                for detail_inst in inst["insts"]: 
                    detail_inst["begin_time"] += begin_time
                    detail_inst["end_time"] += begin_time
                if self.result_json["runtime"] < end_time:
                    self.result_json["runtime"] = end_time
            if i == 0:
                for gate_layer_idx in list_gate_layer_idx:    
                    inst = self.result_json['instructions'][gate_layer_idx]
                    begin_time = self.get_begin_time(gate_layer_idx, inst["dependency"])
                    if inst["type"] == "rydberg":
                        end_time = begin_time + self.architecture.time_rydberg
                    else:
                        end_time = begin_time + (self.architecture.time_1qGate * len(inst["gates"])) + self.common_1q
                    if self.result_json["runtime"] < end_time:
                        self.result_json["runtime"] = end_time
                    inst["begin_time"] = begin_time
                    inst["end_time"] = end_time
            

    def get_begin_time(self, cur_inst_idx: int, dependency: dict):
        begin_time = 0
        for dependency_type in dependency:
            if isinstance(dependency[dependency_type], int):
                inst_idx = dependency[dependency_type]
                if begin_time < self.result_json['instructions'][inst_idx]["end_time"]:
                    begin_time = self.result_json['instructions'][inst_idx]["end_time"]
            else:
                if dependency_type == "site":
                    for inst_idx in dependency[dependency_type]:
                        if self.result_json['instructions'][inst_idx]["type"] == "rearrangeJob":
                            atom_transfer_finish_time = 0
                            for detail_inst in self.result_json['instructions'][inst_idx]["insts"]:
                                inst_type = detail_inst["type"].split(":")[0]
                                if inst_type == "activate":
                                    atom_transfer_finish_time = max(detail_inst["end_time"], atom_transfer_finish_time)
                            atom_transfer_begin_time = 0
                            for detail_inst in self.result_json['instructions'][cur_inst_idx]["insts"]:
                                inst_type = detail_inst["type"].split(":")[0]
                                if inst_type == "deactivate":
                                    atom_transfer_begin_time = max(detail_inst["begin_time"], atom_transfer_begin_time)
                                    break
                            tmp_begin_time = atom_transfer_finish_time - atom_transfer_begin_time
                            if begin_time < tmp_begin_time:
                                begin_time = tmp_begin_time
                        else:
                            if begin_time < self.result_json['instructions'][inst_idx]["end_time"]:
                                begin_time = self.result_json['instructions'][inst_idx]["end_time"]
                else:
                    for inst_idx in dependency[dependency_type]:
                        if begin_time < self.result_json['instructions'][inst_idx]["end_time"]:
                            begin_time = self.result_json['instructions'][inst_idx]["end_time"]
        return begin_time
    
    def get_duration(self, inst: dict):
        list_detail_inst = inst["insts"]
        duration = 0
        for detail_inst in list_detail_inst:
            inst_type = detail_inst["type"].split(":")[0]
            detail_inst["begin_time"] = duration
            if inst_type == "activate" or inst_type == "deactivate":
                duration += self.architecture.time_atom_transfer
                detail_inst["end_time"] = duration
            elif inst_type == "move":
                move_duration = 0
                for row_begin, row_end in zip(detail_inst["row_y_begin"], detail_inst["row_y_end"]):
                    for col_begin, col_end in zip(detail_inst["col_x_begin"], detail_inst["col_x_end"]):
                        tmp = self.architecture.movement_duration(col_begin, row_begin, col_end, row_end)
                        if move_duration < tmp:
                            move_duration = tmp
                detail_inst["end_time"] = move_duration + duration
                duration += move_duration
            else:
                raise ValueError
        return duration

    def _init_aod_heaps(self):
        n_aod = len(self.architecture.dict_AOD)
        self.aod_dependency = [0 for _ in range(n_aod)]
        self.aod_end_time = [(0, i) for i in range(n_aod)]
        if not getattr(self, "zone_affinity", False):
            self.zone_to_heap = None
            self.free_aod_heap = None
            return
        y_zone = []
        for zone in self.architecture.entanglement_zone:
            slm = self.architecture.dict_SLM[zone[0]]
            y_zone.append((slm.location[1], self.architecture.dict_SLM[zone[0]].entanglement_id))
        y_zone = sorted(y_zone, key=lambda x: x[0])
        zone_ids = [z for _, z in y_zone]
        self.zone_to_heap = dict()
        self.free_aod_heap = []
        aod_ids = list(range(n_aod))
        if len(zone_ids) > 0 and len(aod_ids) > 0:
            self.zone_to_heap[zone_ids[0]] = [(0, aod_ids.pop(0))]
        if len(zone_ids) > 1 and len(aod_ids) > 0:
            self.zone_to_heap[zone_ids[1]] = [(0, aod_ids.pop(0))]
        for a in aod_ids:
            heapq.heappush(self.free_aod_heap, (0, a))

    def _pop_aod(self, inst):
        if self.zone_to_heap is None:
            return heapq.heappop(self.aod_end_time)
        target_zone = inst.get("target_zone", None)
        if target_zone is not None and target_zone in self.zone_to_heap and len(self.zone_to_heap[target_zone]) > 0:
            t, aid = heapq.heappop(self.zone_to_heap[target_zone])
            return t, aid
        if self.free_aod_heap and len(self.free_aod_heap) > 0:
            return heapq.heappop(self.free_aod_heap)
        return heapq.heappop(self.aod_end_time)

    def _push_aod(self, aod_id, end_time):
        if self.zone_to_heap is None:
            heapq.heappush(self.aod_end_time, (end_time, aod_id))
            return
        for zid, heap in self.zone_to_heap.items():
            ids = [aid for _, aid in heap]
            if aod_id in ids:
                heapq.heappush(heap, (end_time, aod_id))
                return
        heapq.heappush(self.free_aod_heap, (end_time, aod_id))

    def expand_arrangement(self, inst: dict):
        details = []  # all detailed instructions

        # ---------------------- find out number of cols ----------------------
        all_col_x = [] # all the x coord of qubits
        coords = [] # coords of qubits, shape is same as "begin_locs"
        # these coords are going to be updated as we construct the detail insts

        for locs in inst["begin_locs"]:
            coords_row = []
            for loc in locs:
                exact_location = self.architecture.exact_SLM_location(loc[1], loc[2], loc[3])
                coords_row.append({
                    "id": loc[0],
                    "x": exact_location[0],
                    "y": exact_location[1],
                })

                all_col_x.append(exact_location[0])

            coords.append(coords_row)
        
        init_coords = deepcopy(coords)

        all_col_x = sorted(all_col_x)

        # assign AOD column ids based on all x coords needed
        col_x_to_id = {all_col_x[i]: i for i in range(len(all_col_x))}
        # ---------------------------------------------------------------------
        
        # -------------------- activation and parking -------------------------
        all_col_idx_sofar = [] # which col has been activated
        for row_id, locs in enumerate(inst["begin_locs"]): # each row

            row_y = self.architecture.exact_SLM_location(
                locs[0][1],
                locs[0][2],
                locs[0][3],
            )[1]
            row_loc = [locs[0][1], locs[0][2]]

            shift_back = {
                "type": "move",
                "move_type": "before",
                "row_id": [],
                "row_y_begin": [],
                "row_y_end": [],
                "row_loc_begin": [],
                "row_loc_end": [],
                "col_id": [],
                "col_x_begin": [],
                "col_x_end": [],
                "col_loc_begin": [],
                "col_loc_end": [],
                "begin_coord": deepcopy(coords),
                "end_coord": [],
            }

            activate = {
                "type": "activate",
                "row_id": [row_id, ],
                "row_y": [row_y, ],
                "row_loc": [row_loc, ],
                "col_id": [],
                "col_x": [],
                "col_loc": [],
            }

            for j, loc in enumerate(locs):
                col_x = self.architecture.exact_SLM_location(
                    loc[1],
                    loc[2],
                    loc[3],
                )[0]
                col_loc = [loc[1], loc[3]]
                col_id = col_x_to_id[col_x]
                if col_id not in all_col_idx_sofar:
                    all_col_idx_sofar.append(col_id)
                    activate["col_id"].append(col_id)
                    activate["col_x"].append(col_x)
                    activate["col_loc"].append(col_loc)
                else:
                    shift_back["col_id"].append(col_id)
                    shift_back["col_x_begin"].append(col_x + self.PARKING_DIST)
                    shift_back["col_x_end"].append(col_x)
                    shift_back["col_loc_begin"].append([-1, -1])
                    shift_back["col_loc_end"].append(col_loc)
                    coords[row_id][j]["x"] = col_x
            
            shift_back["end_coord"] = deepcopy(coords)

            if len(shift_back["col_id"]) != 0:
                details.append(shift_back)
            details.append(activate)

            if row_id < len(inst["begin_locs"]) - 1:
                parking = {
                    "type": "move",
                    "move_type": "after",
                    "row_id": [row_id, ],
                    "row_y_begin": [row_y, ],
                    "row_y_end": [row_y + self.PARKING_DIST],
                    "row_loc_begin": [row_loc],
                    "row_loc_end": [[-1, -1]],
                    "col_id": [],
                    "col_x_begin": [],
                    "col_x_end": [],
                    "col_loc_begin": [],
                    "col_loc_end": [],
                    "begin_coord": deepcopy(coords),
                    "end_coord": [],
                }
                for j, loc in enumerate(locs):
                    col_x = self.architecture.exact_SLM_location(
                        loc[1],
                        loc[2],
                        loc[3],
                    )[0]
                    col_loc = [loc[1], loc[3]]
                    col_id = col_x_to_id[col_x]
                    parking["col_id"].append(col_id)
                    parking["col_x_begin"].append(col_x)
                    parking["col_x_end"].append(col_x + self.PARKING_DIST)
                    parking["col_loc_begin"].append(col_loc)
                    parking["col_loc_end"].append([-1, -1])
                    coords[row_id][j]["x"] = parking["col_x_end"][-1]
                    coords[row_id][j]["y"] = parking["row_y_end"][0]
                parking["end_coord"] = deepcopy(coords)
                details.append(parking)
        # ---------------------------------------------------------------------

        # ------------------------- big move ----------------------------------
        big_move = {
            "type": "move:big",
            "move_type": "big",
            "row_id": [],
            "row_y_begin": [],
            "row_y_end": [],
            "row_loc_begin": [],
            "row_loc_end": [],
            "col_id": [],
            "col_x_begin": [],
            "col_x_end": [],
            "col_loc_begin": [],
            "col_loc_end": [],
            "begin_coord": deepcopy(coords),
            "end_coord": [],
        }

        for row_id, (begin_locs, end_locs) in enumerate(zip(
            inst["begin_locs"], inst["end_locs"],
            )):

            big_move["row_id"].append(row_id)
            big_move["row_y_begin"].append(
                coords[row_id][0]["y"]
            )
            if init_coords[row_id][0]["y"] == coords[row_id][0]["y"]:
                big_move["row_loc_begin"].append([begin_locs[0][1], begin_locs[0][2]])
            else:
                big_move["row_loc_begin"].append([-1, -1])

            big_move["row_y_end"].append(
                self.architecture.exact_SLM_location(
                    end_locs[0][1],
                    end_locs[0][2],
                    end_locs[0][3],
                )[1]
            )
            big_move["row_loc_end"].append([end_locs[0][1], end_locs[0][2]])

            for j, (begin_loc, end_loc) in enumerate(zip(begin_locs, end_locs)):
                col_x = self.architecture.exact_SLM_location(
                            begin_loc[1],
                            begin_loc[2],
                            begin_loc[3],
                        )[0]
                col_id = col_x_to_id[col_x]

                if col_id not in big_move["col_id"]:
                    big_move["col_id"].append(col_id)
                    big_move["col_x_begin"].append(coords[row_id][j]["x"])
                    if init_coords[row_id][j]["x"] == coords[row_id][j]["x"]:
                        big_move["col_loc_begin"].append([begin_loc[1], begin_loc[3]])
                    else:
                        big_move["col_loc_begin"].append([-1, -1])
                    big_move["col_x_end"].append(
                        self.architecture.exact_SLM_location(
                            end_loc[1],
                            end_loc[2],
                            end_loc[3],
                        )[0]
                    )
                    big_move["col_loc_end"].append([end_loc[1], end_loc[3]])

                coords[row_id][j]["x"] = self.architecture.exact_SLM_location(
                                            end_loc[1],
                                            end_loc[2],
                                            end_loc[3],
                                        )[0]
                coords[row_id][j]["y"] = self.architecture.exact_SLM_location(
                    end_locs[0][1],
                    end_locs[0][2],
                    end_locs[0][3],
                )[1]
        
        big_move["end_coord"] = deepcopy(coords)
        details.append(big_move)
        # ---------------------------------------------------------------------

        # --------------------------- deactivation ----------------------------
        details.append({
                "type": "deactivate",
                "row_id": [i for i in range(len(inst["begin_locs"]))],
                "col_id": [i for i in range(len(all_col_x))],
            })
        # ---------------------------------------------------------------------
        
        for inst_counter, detail_inst in enumerate(details):
            detail_inst["id"] = inst_counter

        return details
