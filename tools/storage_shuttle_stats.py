#!/usr/bin/env python3
"""
Count storage<->entanglement shuttles per qubit from a compiled code JSON.

Usage:
  python tools/storage_shuttle_stats.py path/to/code.json
"""
import json
import os
import sys


def load_arch(path):
    with open(path, "r") as f:
        spec = json.load(f)
    entanglement_slm = set()
    if "entanglement_zones" in spec:
        for zone in spec["entanglement_zones"]:
            for slm in zone["slms"]:
                entanglement_slm.add(slm["id"])
    storage_slm = set()
    if "storage_zones" in spec:
        for zone in spec["storage_zones"]:
            for slm in zone["slms"]:
                storage_slm.add(slm["id"])
    return entanglement_slm, storage_slm


def resolve_arch_path(code_path, arch_path):
    if os.path.isabs(arch_path):
        return arch_path
    base_dir = os.path.dirname(os.path.dirname(code_path))
    candidate = os.path.join(base_dir, arch_path)
    if os.path.exists(candidate):
        return candidate
    return arch_path


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    code_path = sys.argv[1]
    code = json.load(open(code_path))
    arch_path = resolve_arch_path(code_path, code.get("architecture_spec_path", ""))
    ent_slm, stor_slm = load_arch(arch_path)

    insts = code["instructions"]
    if not insts or insts[0]["type"] != "init":
        raise ValueError("init instruction missing")
    n_q = len(insts[0]["init_locs"])
    locs = [None for _ in range(n_q)]
    for qid, a, r, c in insts[0]["init_locs"]:
        locs[qid] = (a, r, c)

    def zone_of(aid):
        if aid in ent_slm:
            return "ent"
        if aid in stor_slm:
            return "stor"
        return "unknown"

    shuttles = [0 for _ in range(n_q)]
    for inst in insts[1:]:
        if inst["type"] != "rearrangeJob":
            continue
        for qid, a, r, c in inst["end_locs"]:
            prev_zone = zone_of(locs[qid][0])
            new_zone = zone_of(a)
            if prev_zone != new_zone:
                shuttles[qid] += 1
            locs[qid] = (a, r, c)

    total = sum(shuttles)
    movers = [(i, s) for i, s in enumerate(shuttles) if s > 0]
    print(f"Total shuttles stor<->ent: {total}")
    print(f"Qubits with moves: {len(movers)} / {n_q}")
    if movers:
        movers.sort(key=lambda x: x[1], reverse=True)
        print("Top movers (qid: shuttles):", movers[:10])


if __name__ == "__main__":
    main()