# LazySets.jl

[![Build Status](https://github.com/JuliaReach/LazySets.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/JuliaReach/LazySets.jl/actions/workflows/ci.yml?query=branch%3Amaster)
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

This library has been applied in a number of scientic works.

<details>
<summary>List</summary>
  
The articles appear in reverse chronological order.
  
[18] **Reachability of weakly nonlinear systems using Carleman linearization.** Marcelo Forets and Christian Schilling. arXiv preprint arXiv:2108.10390 (2021). Accepted in: 15th International Conference on Reachability Problems (2021). October 25-27 2021, Liverpool, UK.

[17] **Combined Exact and Heuristics Based Approach to Hamiltonian Path Problem Optimization for Route Planning.** Fernando Hernandez, Rafael Sotelo and Marcelo Forets (2021). In Winkenbach, M., Parks, S., & Noszek, J. (Eds.), [Technical Proceedings of the 2021 Amazon Last Mile Routing Research Challenge](https://hdl.handle.net/1721.1/131235) (pp. XXI.1–XXI.12). MIT Libraries.

[16] **ARCH-COMP21 Category Report: Artificial Intelligence and Neural Network Controlled Systems for Continuous and Hybrid Systems Plants.** Taylor T. Johnson, Diego Manzanas Lopez, Luis Benet, Marcelo Forets, Christian Schilling, Radoslav Ivanov, Taylor Carpenter, James Weimer, and Insup Lee. (2021) ARCH21. 8th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH21).

[15] **ARCH-COMP21 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Erika Abraham, Marcelo Forets, Goran Frehse, Daniel Freire, Christian Schilling, Stefan Schupp, and Mark Wetzlinger. (2021) ARCH21. 8th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH21).

[14] **ARCH-COMP21 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Luca Geretti, Julien Alexandre dit Sandretto, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Pieter Collins, Parasara Sridhar Duggirala, Marcelo Forets, Edward Kim, Uziel Linares, David P. Sanders, Christian Schilling, and Mark Wetzlinger. (2021) ARCH21. 8th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH21).
  
[13] **Combining Set Propagation with Finite Element Methods for Time Integration in Transient Solid Mechanics Problems.** Forets, Marcelo, Daniel Freire Caporale, and Jorge M. Pérez Zerpa. arXiv preprint [arXiv:2105.05841](https://arxiv.org/abs/2105.05841) (2021).

[12] **Synthesis of hybrid automata with affine dynamics from time-series data.** García Soto, Miriam, Thomas A. Henzinger, and Christian Schilling. [HSCC](https://hscc.acm.org/2021/) (2021). doi: [10.1145/3447928.3456704](https://doi.org/10.1145/3447928.3456704). [Get pdf from arXiv: 2102.12734](https://arxiv.org/abs/2102.12734)

[11] **Efficient reachability analysis of parametric linear hybrid systems with time-triggered transitions.** Marcelo Forets, Daniel Freire, Christian Schilling, 2020. [arXiv: 2006.12325](https://arxiv.org/abs/2006.12325). [18th ACM-IEEE International Conference on Formal Methods and Models for System Design
](https://ieeexplore.ieee.org/document/9314994). See [conference page](https://iitjammu.ac.in/conferences/memocode2020/index.html).

[10] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Zongnan Bao, Marcelo Forets, Daniel Freire, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling, Stefan Schupp, and Mark Wetzlinger (2020) ARCH20. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH20), vol 74, pages 16--48. [10.29007/7dt2](https://easychair.org/publications/paper/DRpS).

[9] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Luca Geretti, Julien Alexandre dit Sandretto, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Pieter Collins, Marcelo Forets, Daniel Freire, Fabian Immler, Niklas Kochdumper, David P. Sanders and Christian Schilling (2020) ARCH20. 7th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH20), vol 74, pages 49--75. [10.29007/zkf6](https://easychair.org/publications/paper/nrdD).

[8] **Case Study: Reachability Analysis of a unified Combat-Command-and-Control Model.** Sergiy Bogomolov, Marcelo Forets, Kostiantyn Potomkin. *International Conference on Reachability Problems (2020). Lecture Notes in Computer Science, vol 12448.* (2020) doi: [10.1007/978-3-030-61739-4_4](https://dx.doi.org/10.1007/978-3-030-61739-4_4). Presented in the [14th International Conference on Reachability Problems 2020](https://www.irif.fr/~rp2020/). [article](https://link.springer.com/chapter/10.1007/978-3-030-61739-4_4)

[7] **Reachability analysis of linear hybrid systems via block decomposition.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. *IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, 39:11 (2020).* doi: [10.1109/TCAD.2020.3012859](https://dx.doi.org/10.1109/TCAD.2020.3012859). Presented in [Embedded Systems Week 2020](http://esweek.hosting2.acm.org/). [Get pdf from arXiv: 1905.02458](https://arxiv.org/abs/1905.02458).

[6] **Algorithms for verifying deep neural networks.** Liu, C., Arnon, T., Lazarus, C., Barrett, C., & Kochenderfer, M. J. (2019). [Get pdf from arXiv: 1903.06758](https://arxiv.org/abs/1903.06758).

[5] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Rajarshi Ray, Christian Schilling and Stefan Schupp (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 14--40 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP).

[4] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Fabian Immler, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Marcelo Forets, Luca Geretti, Niklas Kochdumper, David P. Sanders and Christian Schilling (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 41--61 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP).

[3] **JuliaReach: a Toolbox for Set-Based Reachability.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. [HSCC](http://hscc2019.eecs.umich.edu/): 22nd ACM International Conference on Hybrid Systems: Computation and Control, see [ACM link here](https://dl.acm.org/citation.cfm?id=3311804). [Get pdf from arXiv: 1901.10736](https://arxiv.org/abs/1901.10736).

[2] **ARCH-COMP18 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Xin Chen, Chuchu Fan, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling and Stefan Schupp (2018) ARCH18. 5th International Workshop on Applied Verification of Continuous and Hybrid Systems, 54: 23–52. doi: [10.29007/73mb](https://dx.doi.org/10.29007/73mb).

[1] **Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Frédéric Viry, Andreas Podelski and Christian Schilling (2018) [HSCC](https://www.hscc2018.deib.polimi.it/) Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control: 41–50. See the [ACM Digital Library link](http://dx.doi.org/10.1145/3178126.3178128), or the [arXiv: 1801.09526](https://arxiv.org/abs/1801.09526).

</details>

## :earth_africa: Ecosystem

`LazySets.jl` is applied in a number of projects in the Julia technical computing stack:

- [ReachabilityAnalysis.jl](https://github.com/JuliaReach/ReachabilityAnalysis.jl) -- Methods to compute sets of states reachable by dynamical systems.
- [NeuralNetworkAnalysis.jl](https://github.com/JuliaReach/NeuralNetworkAnalysis.jl) -- Methods to verify neural network controlled systems using reachability analysis
- [NeuralVerification.jl](https://github.com/sisl/NeuralVerification.jl) -- Methods to soundly verify deep neural networks.
- [InvariantSets.jl](https://github.com/ueliwechsler/InvariantSets.jl) -- Compute, approximate and display invariant sets.
- [OpticSim.jl](https://github.com/microsoft/OpticSim.jl) -- Ray tracing for procedurally generated systems.
- [Photometry.jl](https://github.com/JuliaAstro/Photometry.jl) -- Utilities for characterizing sources in astronomical images.
- [Swalbe.jl](https://github.com/Zitzeronion/Swalbe.jl) -- Simple Julia Lattice Boltzmann Solver for Thin Liquid Films and Droplets

## 👨‍🏫 Workshop at JuliaCon 2021

<details>
<summary>Abstract</summary>

We present [JuliaReach](https://github.com/JuliaReach), a Julia ecosystem to perform reachability analysis of dynamical systems. JuliaReach builds on sound scientific approaches and was, in two occasions (2018 and 2020) the winner of the annual friendly competition on Applied Verification for Continuous and Hybrid Systems ([ARCH-COMP](https://cps-vo.org/group/ARCH)).

The workshop consists of three parts (respectively packages) in [JuliaReach](https://github.com/JuliaReach): our core package for set representations, our main package for reachability analysis, and a new package applying reachability analysis with potential use in domain of control, robotics and autonomous systems.

In the first part we present [LazySets.jl](https://github.com/JuliaReach/LazySets.jl), which provides ways to symbolically represent sets of points as geometric shapes, with a special focus on convex sets and polyhedral approximations. [LazySets.jl](https://github.com/JuliaReach/LazySets.jl) provides methods to apply common set operations, convert between different set representations, and efficiently compute with sets in high dimensions.

In the second part we present [ReachabilityAnalysis.jl](https://github.com/JuliaReach/ReachabilityAnalysis.jl), which provides tools to approximate the set of reachable states of systems with both continuous and mixed discrete-continuous dynamics, also known as hybrid systems. It implements conservative discretization and set-propagation techniques at the state-of-the-art.

In the third part we present [NeuralNetworkAnalysis.jl](https://github.com/JuliaReach/NeuralNetworkAnalysis.jl), which is an application of [ReachabilityAnalysis.jl](https://github.com/JuliaReach/ReachabilityAnalysis.jl) to analyze dynamical systems that are controlled by neural networks. This package can be used to validate or invalidate specifications, for instance about the safety of such systems.

Workshop materials are available here: https://github.com/JuliaReach/JuliaCon-2021-Workshop-Its-All-Set
</details>

[![JuliaCon 2021 video](https://img.youtube.com/vi/P4I7pTvQ4nk/0.jpg)](https://youtu.be/P4I7pTvQ4nk)
