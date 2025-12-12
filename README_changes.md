# ZAC Compiler Tweaks

This README summarizes the compiler-side changes we added, where they live in the code, and how to rerun the experiments.

## Feature flags & locations

- `zone_batched_routing`: routes move layers in zone-sized batches (alternating zone order, MIS per batch).  
  - `zac/router/router.py` ~50–200: `_process_movement_batches`, forward routing, reverse handling.

- `interleave_reverse`: immediately drains reverse moves after each 2Q stage using the same batching.  
  - `zac/router/router.py` ~80–160: reverse graph build + stats logging.

- `locality_reorder` and `locality_cluster_split`: reorder 2Q layer gates using previous-stage Rydberg placements; optionally split oversized layers into contiguous clusters.  
  - `zac/scheduler/scheduler.py` ~150–250: sort keys, alternating zone priority, cluster split, TRACE logs.

- `keep_hot`, `keep_hot_penalty`, `keep_hot_horizon`: bias placement to keep reuse qubits in the same entanglement zone across nearby layers (including gaps up to the horizon).  
  - `zac/zac.py` ~35–87: flag parsing; ~246–354: `extend_keep_hot_gap` marks reuse across gaps.  
  - `zac/placer/placer.py` ~98–105: pass flags into placer.  
  - `zac/placer/vmplacer.py` ~8–20, ~269–281: cost penalty when a reuse qubit would switch zones.

- `zone_affinity` (multi-AOD): tie AODs to fixed zones (top/bot; extra AOD free if odd count).  
  - `zac/router/router.py` ~200–260: AOD heap init/pop/push honors affinity.

- Routing assertion removal to allow keep-hot cases.  
  - `zac/router/router.py`: legacy storage/entanglement ID assertion removed near routing entry.

- Config Hooks: experiment JSONs can override flags.  
  - `zac/zac.py` ~35–90: flags; `run.py` merges `policy_file`/`policies` into each experiment.

## Helper tools

- Zone switch counter: `tools/zone_switch_rate.py`  
  `python tools/zone_switch_rate.py result/.../code/ising_n98_transpiled_code.json`

- Storage shuttle counter: `tools/storage_shuttle_stats.py`  
  `python tools/storage_shuttle_stats.py result/.../code/ising_n98_transpiled_code.json`

- QAOA generator: `tools/gen_qaoa_line_100.py` (builds QAOA-line-100 circuits).

## Running experiments

1) Pick an experiment JSON (e.g., `exp_setting/report/ising98_baseline.json`).  
2) Set policies in `exp_setting/report/ising98_policies.json` (e.g., toggle `keep_hot`, `locality_reorder`, etc.), or point the exp JSON to a different policy file.  
3) Run:
```
python run.py exp_setting/report/ising98_baseline.json
```
Results go under `result/.../` (code/fidelity/time JSONs).


## What changed vs. baseline

- Added optional zone-batched routing + reverse interleave to reduce AOD contention/idle.
- Added locality-aware gate reorder/cluster split to cluster 2Q gates by prior placements.
- Added keep-hot bias with horizon to cut storage round-trips and zone hops for reuse qubits.
- Added zone affinity for multi-AOD setups.
- Added diagnostics/tools for zone switches and shuttles.

All flags default to baseline behavior when off. Physical parameters such as durations/fidelities/site counts untouched.
