# Nanoscale Channel Flow & Reaction Simulator

A minimal Octave + C++ toolbox for computing fluid flow in narrow channels and the steady‑state concentration of a reacting species. Some scripts may not work on on some computers as this was written to be compatible with macOS Seqouia 15.4.1 and has shown many bugs given any update. Basically this project is a case study in why apple should be banned from all scientific research.

---

## Requirements

* Octave 9.12 - Simulation code
* Octave 9.12 or greater - Plotting code
* A C++17‑compatible compiler configured for `oct`

---

## Quick Start

```bash
# clone this repository
$ git clone https://github.com/USER/REPO.git
$ cd REPO

# compile the C++ sources (Unix/macOS)
$ make -C src/cppfiles            # creates *.mex* binaries next to sources

# open Octave and run:
>> cd('src/cppfiles');
>> make   % repeat for other *.cpp files

# add everything to MATLAB path and launch the example sweep
>> addpath(genpath(pwd));
>> main
```

Results (fields, metrics, plots) are written to `output/` and `figures/`.

---

## Folder Layout

```
Helpers_analytical/   # Closed‑form solutions for validation
Helpers_Plotting/     # Scripts to load output and create figures/metrics
src/
  ├─ cppfiles/        # Finite‑difference & LBM kernels (C++, Oct)
  └─ mfiles/          # MATLAB solvers, wrappers, utilities
GeometryGenerator.m   # Build 2‑D channel geometry
LBMSimulator.m        # Flow solver front‑end
ReactionDiffusionSimulator.m  # Advection–diffusion–reaction solver
main.m                # Example batch run over Pe & Da
Plotting_finals.m     # Post‑processing entry point
```

---

## Credits

Some functions inside `src/cppfiles/` and `src/mfiles/` were from code written by **Professor Jesper Schmidt Hansen**. His contributions are gratefully acknowledged.

---

## License

This project is released under the MIT License (see `LICENSE`).
