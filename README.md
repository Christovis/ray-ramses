# Ray-RAMSES

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-2007.03042%20-green.svg)](https://arxiv.org/abs/1601.02012)

A open-source ray tracing code to compute integrated cosmological observables on the fly in AMR N-body simulations. Unlike conventional ray tracing techniques, our code takes full advantage of the time and spatial resolution attained by the N-body simulation by computing the integrals along the line of sight on a cell-by-cell basis through the AMR simulation grid. Moroever, since it runs on the fly in the N-body run, our code can produce maps of the desired observables without storing large (or any) amounts of data for post-processing. The ray tracing methodology presented here can be used in several cosmological analysis such as Sunyaev-Zel’dovich and integrated Sachs-Wolfe effect studies as well as modified gravity. Our code can also be used in cross-checks of the more conventional methods, which can be important in tests of theory systematics in preparation for upcoming large scale structure surveys.

### Build/Install
The code has ben tested and run on [COSMA](https://www.dur.ac.uk/icc/cosma/), which uses CentOS 7.6 with a 3.10 Linux kernel.

**Prerequisites:**
* A Fortran compiler
* CMake
* Intel MPI (2018)
* IntelComp (2018)

To build the Fortran BMI bindings from source with cmake, run

```
$ cd ./src/bin
$ make
```

You can quickly test your installation by executing:
```
$ cd bin
$ make
$ cd ..
$ bin/ramses1d namelist/tube1d.nml
```

Initial condition can be generated with [2LPTic](https://arxiv.org/abs/astro-ph/0606505).

The simulation results can be analysed through [Astrild](https://github.com/Christovis/astrild).
