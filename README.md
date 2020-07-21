# LazySets.jl

[![Build Status](https://travis-ci.org/JuliaReach/LazySets.jl.svg?branch=master)](https://travis-ci.org/JuliaReach/LazySets.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliareach.github.io/LazySets.jl/dev/)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/JuliaReach/LazySets.jl/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/105701832.svg)](https://zenodo.org/badge/latestdoi/105701832)
[![Code coverage](http://codecov.io/github/JuliaReach/LazySets.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaReach/LazySets.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

`LazySets` is a [Julia](http://julialang.org) package for calculus with convex sets.

## Resources

- [Manual](http://juliareach.github.io/LazySets.jl/dev/)
- [Contributing](https://juliareach.github.io/LazySets.jl/dev/about/#Contributing-1)
- [Release notes of tagged versions](https://github.com/JuliaReach/LazySets.jl/releases)
- [Release notes of the development version](https://github.com/JuliaReach/LazySets.jl/wiki/Release-log-tracker)
- [Publications](https://juliareach.github.io/Reachability.jl/dev/publications/)
- [Developers](https://juliareach.github.io/LazySets.jl/dev/about/#Credits-1)

## Installing

`LazySets.jl` is a registered Julia package and [as such you can install it](https://julialang.github.io/Pkg.jl/v1/managing-packages/) by activating the `pkg` mode (type `]`, and to leave it, type `<backspace>`),
followed by

```julia
pkg> add LazySets
```

See the [Getting started](https://juliareach.github.io/LazySets.jl/dev/man/getting_started/) section of the manual for other options.

```julia
pkg> add LazySets
```

## Ecosystem

`LazySets.jl` is applied in a number of projects in the Julia technical computing stack:

- [ReachabilityAnalysis.jl](https://github.com/JuliaReach/ReachabilityAnalysis.jl) -- Methods to compute sets of states reachable by dynamical systems
- [NeuralVerification.jl](https://github.com/sisl/NeuralVerification.jl) -- Methods to soundly verify deep neural networks

## Publications

This library has been applied in a number of scientic works. In reverse chronological order,

- ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics. Fabian Immler, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Marcelo Forets, Luca Geretti, Niklas Kochdumper, David P. Sanders and Christian Schilling (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 41–61 doi: 10.29007/bj1w. Packages: Reachability.jl.

- ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics. Matthias Althoff, Stanley Bak, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Rajarshi Ray, Christian Schilling and Stefan Schupp (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 14–40 doi: 10.29007/bj1w. Packages: Reachability.jl.

- Reachability analysis of linear hybrid systems via block decomposition. Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. Get pdf from arXiv: 1905.02458.

- JuliaReach: a Toolbox for Set-Based Reachability. Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. Published in Proceedings of HSCC'19: 22nd ACM International Conference on Hybrid Systems: Computation and Control (HSCC'19), see ACM link here. Get pdf from arXiv: 1901.10736.

- ARCH-COMP18 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics. Matthias Althoff, Stanley Bak, Xin Chen, Chuchu Fan, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling and Stefan Schupp (2018) ARCH18. 5th International Workshop on Applied Verification of Continuous and Hybrid Systems, 54: 23–52. doi: 10.29007/73mb. Packages: Reachability.jl.

- Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices. Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Frédéric Viry, Andreas Podelski and Christian Schilling (2018) HSCC'18 Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control: 41–50. See the ACM Digital Library link, or the arXiv: 1801.09526. Packages: LazySets.jl and Reachability.jl.
