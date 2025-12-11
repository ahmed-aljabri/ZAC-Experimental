# ZAC (fork) — Compiler Experiments for Zoned Neutral-Atom Architectures

> This repository is a **fork** of the original ZAC compiler by **UCLA-VAST**.  
> Upstream: https://github.com/UCLA-VAST/ZAC  
> This fork is not affiliated with or endorsed by the upstream authors.

This fork preserves ZAC’s physics models and core compilation flow, and adds **optional, compiler-side policies** for research experiments (all disabled by default). Examples include zone-batched routing, locality-aware in-layer reordering, reverse interleaving, keep-hot bias with a horizon, and zone-affinity for multi-AOD.

For a detailed list of modifications and how to use them, see **`README_changes.md`**.

---

## License & Attribution

- Licensed under the **BSD-3-Clause** license; see **`LICENSE`** at the repository root.
---


## Getting Started
- Use Python 3 in a virtual environment.
- Install dependencies: `python3 -m pip install -r requirements.txt`
- (Optional) Install `ffmpeg` if you want animation generation.

## Running the Compiler
Run an experiment JSON (e.g., Ising98 baseline):
```
python run.py exp_setting/report/ising98_baseline.json
```
Results (code/fidelity/time) are written under `result/...`.

To toggle the experimental policies, edit the policy file referenced by the experiment (e.g., `exp_setting/report/ising98_policies.json`) or supply a different policy file. Flags include zone-batched routing, locality-aware reorder, reverse interleave, keep-hot with horizon, and zone affinity. Details are in `README_changes.md`.

## Repository Structure (high level)
- `zac/` – compiler source (scheduler, placer, router, verifier, animator).
- `exp_setting/` – experiment and policy JSONs.
- `hardware_spec/` – architecture specs.
- `benchmark/` – example circuits.
- `result/` – outputs (code, fidelity, time, animations).
- `tools/` – diagnostics (zone switches, shuttle counts, QAOA generator, etc.).

## References
Original ZAC publication:
```
@inproceedings{lin2025reuse,
  title={Reuse-aware compilation for zoned quantum architectures based on neutral atoms},
  author={Lin, Wan-Hsuan and Tan, Daniel Bochen and Cong, Jason},
  booktitle={2025 IEEE International Symposium on High Performance Computer Architecture (HPCA)},
  pages={127--142},
  year={2025},
  organization={IEEE}
}
```
