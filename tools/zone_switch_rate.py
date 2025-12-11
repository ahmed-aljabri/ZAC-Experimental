#!/usr/bin/env python3
"""
Compute how often qubits switch entanglement zones between consecutive Rydberg layers.

Usage:
  python tools/zone_switch_rate.py path/to/code.json
"""
import json
import os
import sys


def load_arch(path):
    with open(path, "r") as f:
        spec = json.load(f)
    slm_to_zone = {}
    if "storage_zones" in spec:
        for zone in spec["storage_zones"]:
            for slm in zone["slms"]:
                slm_to_zone[slm["id"]] = -1
    if "entanglement_zones" in spec:
        for zone in spec["entanglement_zones"]:
            zid = zone["zone_id"]
            for slm in zone["slms"]:
                slm_to_zone[slm["id"]] = zid
    return slm_to_zone


def resolve_arch_path(code_path, arch_path):
    if os.path.isabs(arch_path):
        return arch_path
    base_dir = os.path.dirname(os.path.dirname(code_path))
    candidate = os.path.join(base_dir, arch_path)
    if os.path.exists(candidate):
        return candidate
    # fall back to cwd
    return arch_path


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    code_path = sys.argv[1]
    with open(code_path, "r") as f:
        code = json.load(f)
    arch_path = resolve_arch_path(code_path, code.get("architecture_spec_path", ""))
    slm_to_zone = load_arch(arch_path)

    insts = code["instructions"]
    if not insts or insts[0]["type"] != "init":
        raise ValueError("First instruction must be init")
    n_qubits = len(insts[0]["init_locs"])
    locs = [None for _ in range(n_qubits)]
    for qid, a, r, c in insts[0]["init_locs"]:
        locs[qid] = (a, r, c)

    last_zone = [None for _ in range(n_qubits)]
    switch_count = [0 for _ in range(n_qubits)]
    events = [0 for _ in range(n_qubits)]

    for inst in insts[1:]:
        if inst["type"] == "rearrangeJob":
            for qid, a, r, c in inst["end_locs"]:
                locs[qid] = (a, r, c)
        elif inst["type"] == "rydberg":
            seen_in_this_inst = set()
            for gate in inst["gates"]:
                for qid in (gate["q0"], gate["q1"]):
                    if qid in seen_in_this_inst:
                        continue
                    seen_in_this_inst.add(qid)
                    a = locs[qid][0]
                    zone = slm_to_zone.get(a, None)
                    if zone is None or zone == -1:
                        # unexpected placement
                        continue
                    if last_zone[qid] is not None:
                        events[qid] += 1
                        if zone != last_zone[qid]:
                            switch_count[qid] += 1
                    else:
                        # first time this qubit participates
                        pass
                    last_zone[qid] = zone

    total_events = sum(events)
    total_switches = sum(switch_count)
    rate = (total_switches / total_events) if total_events else 0.0

    print(f"Total switch events: {total_switches}/{total_events} (rate={rate:.4f})")
    print("Per-qubit switches (qid: switches/events):")
    for qid, (sw, ev) in enumerate(zip(switch_count, events)):
        if ev:
            print(f"  {qid}: {sw}/{ev} ({sw/ev:.4f})")


if __name__ == "__main__":
    main()