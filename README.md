# Nanoscale Channel Flow & Reaction Simulator

A minimal MATLAB + C++ MEX toolbox for computing fluid flow in narrow channels and the steady‑state concentration of a reacting species.

---

## Requirements

* MATLAB R2024a (or newer)
* A C++17‑compatible compiler configured for `mex`

---

## Quick Start

```bash
# clone this repository
$ git clone https://github.com/USER/REPO.git
$ cd REPO

# compile the C++ sources (Unix/macOS)
$ make -C src/cppfiles            # creates *.mex* binaries next to sources

# on Windows, open MATLAB and run:
>> cd('src/cppfiles');
>> mex -R2018a -largeArrayDims computeQMap_mex.cpp   % repeat for other *.cpp files

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
  ├─ cppfiles/        # Finite‑difference & LBM kernels (C++, MEX)
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
