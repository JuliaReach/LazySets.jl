# LazySets.jl

[![Build Status](https://travis-ci.org/JuliaReach/LazySets.jl.svg?branch=master)](https://travis-ci.org/JuliaReach/LazySets.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliareach.github.io/LazySets.jl/dev/)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/JuliaReach/LazySets.jl/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/105701832.svg)](https://zenodo.org/badge/latestdoi/105701832)
[![Code coverage](http://codecov.io/github/JuliaReach/LazySets.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaReach/LazySets.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

`LazySets` is a [Julia](http://julialang.org) package for calculus with convex sets.

## 🎯 Resources

- [Manual](http://juliareach.github.io/LazySets.jl/dev/)
- [Contributing](https://juliareach.github.io/LazySets.jl/dev/about/#Contributing-1)
- [Release notes of tagged versions](https://github.com/JuliaReach/LazySets.jl/releases)
- [Release notes of the development version](https://github.com/JuliaReach/LazySets.jl/wiki/Release-log-tracker)
- [Developers](https://juliareach.github.io/LazySets.jl/dev/about/#Credits-1)

## 💾 Installing

`LazySets.jl` is a registered Julia package and [as such you can install it](https://julialang.github.io/Pkg.jl/v1/managing-packages/) by activating the `pkg` mode (type `]`, and to leave it, type `<backspace>`),
followed by

```julia
pkg> add LazySets
```

See the [Getting started](https://juliareach.github.io/LazySets.jl/dev/man/getting_started/) section of the manual for other options.

## :blue_book: Publications

This library has been applied in a number of scientic works. In reverse chronological order,

[11] **Efficient reachability analysis of parametric linear hybrid systems with time-triggered transitions.** Marcelo Forets, Daniel Freire, Christian Schilling, 2020. [arXiv: 2006.12325](https://arxiv.org/abs/2006.12325). Accepted in [18th ACM-IEEE International Conference on Formal Methods and Models for System Design
](https://iitjammu.ac.in/conferences/memocode2020/index.html).

[10] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Zongnan Bao, Marcelo Forets, Daniel Freire, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling, Stefan Schupp, and Mark Wetzlinger (2020) ARCH20. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH20), vol 74, pages 16--48. [10.29007/7dt2](https://easychair.org/publications/paper/DRpS).

[9] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Luca Geretti, Julien Alexandre dit Sandretto, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Pieter Collins, Marcelo Forets, Daniel Freire, Fabian Immler, Niklas Kochdumper, David P. Sanders and Christian
Schilling (2020) ARCH20. To appear in 7th International Workshop on Applied Verification of Continuous and Hybrid Systems. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH20), vol 74, pages 49--75. [10.29007/zkf6](https://easychair.org/publications/paper/nrdD).

[8] **Case Study: Reachability Analysis of a unified Combat-Command-and-Control Model.** Sergiy Bogomolov, Marcelo Forets, Kostiantyn Potomkin. Accepted in [14th International Conference on Reachability Problems 2020](https://www.irif.fr/~rp2020/). [article](https://link.springer.com/chapter/10.1007/978-3-030-61739-4_4)

[7] **Reachability analysis of linear hybrid systems via block decomposition.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. [Get pdf from arXiv: 1905.02458](https://arxiv.org/abs/1905.02458). Accepted in [Embedded Systems Week 2020](http://esweek.hosting2.acm.org/).

[6] **Algorithms for verifying deep neural networks.** Liu, C., Arnon, T., Lazarus, C., Barrett, C., & Kochenderfer, M. J. (2019). [Get pdf from arXiv: 1903.06758](https://arxiv.org/abs/1903.06758).

[5] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Rajarshi Ray, Christian Schilling and Stefan Schupp (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 14--40 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP).

[4] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Fabian Immler, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Marcelo Forets, Luca Geretti, Niklas Kochdumper, David P. Sanders and Christian Schilling (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 41--61 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP).

[3] **JuliaReach: a Toolbox for Set-Based Reachability.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. Published in Proceedings of [HSCC'19](http://hscc2019.eecs.umich.edu/): 22nd ACM International Conference on Hybrid Systems: Computation and Control (HSCC'19), see [ACM link here](https://dl.acm.org/citation.cfm?id=3311804). [Get pdf from arXiv: 1901.10736](https://arxiv.org/abs/1901.10736).

[2] **ARCH-COMP18 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Xin Chen, Chuchu Fan, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling and Stefan Schupp (2018) ARCH18. 5th International Workshop on Applied Verification of Continuous and Hybrid Systems, 54: 23–52. doi: 10.29007/73mb.

[1] **Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Frédéric Viry, Andreas Podelski and Christian Schilling (2018) [HSCC'18](https://www.hscc2018.deib.polimi.it/) Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control: 41–50. See the [ACM Digital Library link](http://dx.doi.org/10.1145/3178126.3178128), or the [arXiv: 1801.09526](https://arxiv.org/abs/1801.09526).

## :earth_africa: Ecosystem

`LazySets.jl` is applied in a number of projects in the Julia technical computing stack:

- [ReachabilityAnalysis.jl](https://github.com/JuliaReach/ReachabilityAnalysis.jl) -- Methods to compute sets of states reachable by dynamical systems.
- [NeuralVerification.jl](https://github.com/sisl/NeuralVerification.jl) -- Methods to soundly verify deep neural networks.
- [Photometry.jl](https://github.com/JuliaAstro/Photometry.jl) -- Utilities for characterizing sources in astronomical images.
- [InvariantSets.jl](https://github.com/ueliwechsler/InvariantSets.jl) -- Compute, approximate and display invariant sets.
