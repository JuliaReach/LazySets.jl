# LazySets.jl
## *Scalable Symbolic-Numeric Set Computations*

| **Introduction & Documentation** |**Status** |**Community** |**Version-specific Citation** | **License** |
|:-----------------:|:---------------:|:------------:|:----------------------------:|:-----------:|
| [![paper][paper-img]][paper-url] [![docs-dev][dev-img]][dev-url] | [![CI][ci-img]][ci-url] [![codecov][cov-img]][cov-url] | [![zulip][chat-img]][chat-url] | [![zenodo][doi-img]][doi-url] | [![license][lic-img]][lic-url] |

[paper-img]: https://proceedings.juliacon.org/papers/10.21105/jcon.00097/status.svg
[paper-url]: https://doi.org/10.21105/jcon.00097
[dev-img]: https://img.shields.io/badge/docs-latest-blue.svg
[dev-url]: https://juliareach.github.io/LazySets.jl/dev/
[ci-img]: https://github.com/JuliaReach/LazySets.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/JuliaReach/LazySets.jl/actions/workflows/ci.yml
[cov-img]: https://codecov.io/github/JuliaReach/LazySets.jl/coverage.svg
[cov-url]: https://app.codecov.io/github/JuliaReach/LazySets.jl
[chat-img]: https://img.shields.io/badge/zulip-join_chat-brightgreen.svg
[chat-url]: https://julialang.zulipchat.com/#narrow/stream/278609-juliareach
[doi-img]: https://zenodo.org/badge/105701832.svg
[doi-url]: https://zenodo.org/badge/latestdoi/105701832
[lic-img]: https://img.shields.io/github/license/mashape/apistatus.svg
[lic-url]: https://github.com/JuliaReach/LazySets.jl/blob/master/LICENSE

## ❓ Introduction

The following article showcases the basic functionality, highlighting some of the key design choices:

> Forets, Marcelo, and Christian Schilling. *LazySets.jl: Scalable Symbolic-Numeric Set Computations.* [Proceedings of the JuliaCon Conferences](https://doi.org/10.21105/jcon.00097) (2021).

See [below](#-how-to-cite) for how to cite it.

## 🎯 Resources

- [Manual](https://juliareach.github.io/LazySets.jl/dev/)
- [Release notes of the development version](https://github.com/JuliaReach/LazySets.jl/wiki/Release-log-tracker)
- [Release notes of previous versions](https://github.com/JuliaReach/LazySets.jl/releases)
- [How to contribute](https://juliareach.github.io/LazySets.jl/dev/about/#Contributing-1)
- [Developers team](https://juliareach.github.io/LazySets.jl/dev/about/#Credits-1)

## 💾 Installing

`LazySets.jl` is a registered Julia package and [as such you can install it](https://julialang.github.io/Pkg.jl/v1/managing-packages/) by activating the `pkg` mode (type `]`, and to leave it, type `<backspace>`),
followed by

```julia
pkg> add LazySets
```

See the [Getting started](https://juliareach.github.io/LazySets.jl/dev/man/getting_started/) section of the manual for other options.

## 📘 Publications

This library has been applied in a number of scientific works.

<details>
<summary>Click to see the full list of publications that use LazySets.</summary>

The articles appear in reverse chronological order.

[26] **ARCH-COMP22 category report: Continuous and hybrid systems with linear continuous dynamics.** Matthias Althoff, Marcelo Forets, Christian Schilling, and Mark Wetzlinger. (2022). 9th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH22), vol 90, pp. 58-85. [doi: 10.29007/mmzc](https://doi.org/10.29007/mmzc).

[25] **ARCH-COMP22 category report: Continuous and hybrid systems with nonlinear dynamics.** Luca Geretti, Julien Alexandre Dit Sandretto, Matthias Althoff, Luis Benet, Pieter Collins, Parasara Sridhar Duggirala, Marcelo Forets, Edward Kim, Stefan Mitsch, Christian Schilling, and Mark Wetzlinger. (2022). 9th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH22), vol 90, pp. 86-112. [doi: 10.29007/fnzc](https://doi.org/10.29007/fnzc).

[24] **ARCH-COMP22 category report: Artificial intelligence and neural network control systems for continuous and hybrid systems plants.** Diego Manzanas Lopez, Matthias Althoff, Luis Benet, Xin Chen, Jiameng Fan, Marcelo Forets, Chao Huang, Taylor T. Johnson, Tobias Ladner, Wenchao Li, Christian Schilling, and Qi Zhu. (2022). 9th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH22), vol 90, pp. 142-184. [doi: 10.29007/wfgr](https://doi.org/10.29007/wfgr).

[23] **Synthesis of parametric hybrid automata from time series.** Miriam García Soto, Thomas A. Henzinger, and Christian Schilling (2022). Proceedings of the [20th International Symposium on Automated Technology for Verification and Analysis](https://atva-conference.org/2022/), LNCS, vol. 13505, pp. 337-353. [doi: 10.1007/978-3-031-19992-9_22](https://doi.org/10.1007/978-3-031-19992-9_22), [arXiv: 2208.06383](https://arxiv.org/abs/2208.06383).

[22] **Decomposing reach set computations with low-dimensional sets and high-dimensional matrices (extended version).** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Andreas Podelski, and Christian Schilling (2022). [Information and Computation](https://www.sciencedirect.com/journal/information-and-computation), vol. 289. [doi: 10.1016/j.ic.2022.104937](https://doi.org/10.1016/j.ic.2022.104937).

[21] **Conservative Time Discretization: A Comparative Study.** Marcelo Forets and Christian Schilling (2022). Proceedings of the [17th International Conference on integrated Formal Methods (iFM)](https://ifm22.si.usi.ch/), LNCS, vol. 13274, pp. 149-167. [doi: 10.1007/978-3-031-07727-2_9](https://doi.org/10.1007/978-3-031-07727-2_9), [arXiv: 2111.01454](https://arxiv.org/abs/2111.01454).

[20] **Verification of Neural-Network Control Systems by Integrating Taylor Models and Zonotopes.** Christian Schilling, Marcelo Forets, and Sebastián Guadalupe (2022). Proceedings of the [36th Conference on Artificial Intelligence (AAAI)](https://aaai.org/Conferences/AAAI-22/). [doi: 10.1609/aaai.v36i7.20790](https://doi.org/10.1609/aaai.v36i7.20790).

[19] **Combining Set Propagation with Finite Element Methods for Time Integration in Transient Solid Mechanics Problems.** Marcelo Forets, Daniel Freire Caporale, and Jorge M. Pérez Zerpa (2022). [Computers & Structures](https://www.sciencedirect.com/journal/computers-and-structures), vol 259. [doi: 10.1016/j.compstruc.2021.106699](https://doi.org/10.1016/j.compstruc.2021.106699), [arXiv: 2105.05841](https://arxiv.org/abs/2105.05841).

[18] **LazySets.jl: Scalable Symbolic-Numeric Set Computations.** Marcelo Forets and Christian Schilling (2021). [Proceedings of the JuliaCon Conferences](https://proceedings.juliacon.org/). [doi: 10.21105/jcon.00097](https://doi.org/10.21105/jcon.00097).

[17] **Reachability of weakly nonlinear systems using Carleman linearization.** Marcelo Forets and Christian Schilling (2021). Proceedings of the [15th International Conference on Reachability Problems (RP)](https://rp2021.csc.liv.ac.uk/), LNCS, vol. 13035, pp. 85-99. [doi: 10.1007/978-3-030-89716-1_6](https://doi.org/10.1007/978-3-030-89716-1_6), [arXiv: 2108.10390](https://arxiv.org/abs/2108.10390).

[16] **Combined Exact and Heuristics Based Approach to Hamiltonian Path Problem Optimization for Route Planning.** Fernando Hernandez, Rafael Sotelo, and Marcelo Forets (2021). Technical Proceedings of the [2021 Amazon Last Mile Routing Research Challenge](https://hdl.handle.net/1721.1/131235).

[15] **ARCH-COMP21 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Erika Abraham, Marcelo Forets, Goran Frehse, Daniel Freire, Christian Schilling, Stefan Schupp, and Mark Wetzlinger. (2021). 8th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH21), vol 80, pp. 1-31. [doi: 10.29007/lhbw](https://doi.org/10.29007/lhbw).

[14] **ARCH-COMP21 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Luca Geretti, Julien Alexandre dit Sandretto, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Pieter Collins, Parasara Sridhar Duggirala, Marcelo Forets, Edward Kim, Uziel Linares, David P. Sanders, Christian Schilling, and Mark Wetzlinger. (2021). 8th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH21), vol 80, pp. 32-54. [doi: 10.29007/2jw8](https://doi.org/10.29007/2jw8).

[13] **ARCH-COMP21 Category Report: Artificial Intelligence and Neural Network Controlled Systems for Continuous and Hybrid Systems Plants.** Taylor T. Johnson, Diego Manzanas Lopez, Luis Benet, Marcelo Forets, Christian Schilling, Radoslav Ivanov, Taylor Carpenter, James Weimer, and Insup Lee. (2021). 8th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH21), vol 80, pp. 90-119. [doi: 10.29007/kfk9](https://doi.org/10.29007/kfk9).

[12] **Synthesis of hybrid automata with affine dynamics from time-series data.** Miriam García Soto, Thomas A. Henzinger, and Christian Schilling (2021). [24th International Conference on Hybrid Systems: Computation and Control (HSCC)](https://hscc.acm.org/2021/). [doi: 10.1145/3447928.3456704](https://doi.org/10.1145/3447928.3456704), [arXiv: 2102.12734](https://arxiv.org/abs/2102.12734).

[11] **Algorithms for verifying deep neural networks.** Changliu Liu, Tomer Arnon, Christopher Lazarus, Christopher A. Strong, Clark W. Barrett, and Mykel J. Kochenderfer (2021). [Foundations and Trends in Optimization](https://www.nowpublishers.com/OPT), vol 4, pp. 244-404. [doi: 10.1561/2400000035](https://doi.org/10.1561/2400000035), [arXiv: 1903.06758](https://arxiv.org/abs/1903.06758).

[10] **Efficient reachability analysis of parametric linear hybrid systems with time-triggered transitions.** Marcelo Forets, Daniel Freire, and Christian Schilling (2020). Proceedings of the [18th International Conference on Formal Methods and Models for System Design (MEMOCODE)](https://iitjammu.ac.in/conferences/memocode2020/), pp. 137-142. [doi: 10.1109/MEMOCODE51338.2020.9314994](https://doi.org/10.1109/MEMOCODE51338.2020.9314994), [arXiv: 2006.12325](https://arxiv.org/abs/2006.12325).

[9] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Zongnan Bao, Marcelo Forets, Daniel Freire, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling, Stefan Schupp, and Mark Wetzlinger (2020). 7th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH20), vol 74, pp. 16-48. [doi: 10.29007/7dt2](https://doi.org/10.29007/7dt2).

[8] **ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Luca Geretti, Julien Alexandre dit Sandretto, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Pieter Collins, Marcelo Forets, Daniel Freire, Fabian Immler, Niklas Kochdumper, David P. Sanders, and Christian Schilling (2020). 7th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH20), vol 74, pp. 49-75. [doi: 10.29007/zkf6](https://doi.org/10.29007/zkf6).

[7] **Case Study: Reachability Analysis of a unified Combat-Command-and-Control Model.** Sergiy Bogomolov, Marcelo Forets, and Kostiantyn Potomkin (2020). Proceedings of the [14th International Conference on Reachability Problems (RP)](https://www.irif.fr/~rp2020/), LNCS, vol 12448, pp. 52-66. [doi: 10.1007/978-3-030-61739-4_4](https://doi.org/10.1007/978-3-030-61739-4_4).

[6] **Reachability analysis of linear hybrid systems via block decomposition.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, and Christian Schilling (2020). IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, vol. 39, pp. 4018-4029. Presented at [Embedded Systems Week](https://esweek.org/) 2020. [doi: 10.1109/TCAD.2020.3012859](https://doi.org/10.1109/TCAD.2020.3012859), [arXiv: 1905.02458](https://arxiv.org/abs/1905.02458).

[5] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Rajarshi Ray, Christian Schilling, and Stefan Schupp (2019). 6th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH19), vol 61, pp. 14-40. [doi: 10.29007/bj1w](https://doi.org/10.29007/bj1w).

[4] **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Fabian Immler, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Marcelo Forets, Luca Geretti, Niklas Kochdumper, David P. Sanders, and Christian Schilling (2019). 6th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH19), vol 61, pp. 41-61. [doi: 10.29007/m75b](https://doi.org/10.29007/m75b).

[3] **JuliaReach: a Toolbox for Set-Based Reachability.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling (2019). Proceedings of the 22nd International Conference on Hybrid Systems: Computation and Control (HSCC), pp. 39-44. [doi: 10.1145/3302504.3311804](https://doi.org/10.1145/3302504.3311804), [arXiv: 1901.10736](https://arxiv.org/abs/1901.10736).

[2] **ARCH-COMP18 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Xin Chen, Chuchu Fan, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling, and Stefan Schupp (2018). 5th [International Workshop on Applied Verification of Continuous and Hybrid Systems](https://cps-vo.org/group/ARCH/) (ARCH18), vol 54, pp. 23–52. [doi: 10.29007/73mb](https://doi.org/10.29007/73mb).

[1] **Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Frédéric Viry, Andreas Podelski, and Christian Schilling (2018). Proceedings of the [21st International Conference on Hybrid Systems: Computation and Control (HSCC)](https://www.hscc2018.deib.polimi.it/), pp. 41–50. [doi: 10.1145/3178126.3178128](https://doi.org/10.1145/3178126.3178128), [arXiv: 1801.09526](https://arxiv.org/abs/1801.09526).

</details>

## 🗺 Ecosystem

Several projects in the Julia technical computing stack use this library.

<details>
<summary>Click to see the full list of Julia packages that use LazySets.</summary>

- [ConvexBodyProximityQueries.jl](https://github.com/arlk/ConvexBodyProximityQueries.jl) -- Proximity computation between convex bodies in 2D/3D.
- [ClosedLoopReachability.jl](https://github.com/JuliaReach/ClosedLoopReachability.jl) -- Reachability analysis for closed-loop control systems.
- [HySynthParametric](https://github.com/HySynth/HySynthParametric) -- Synthesis of parametric linear hybrid automata.
- [IntervalLinearAlgebra.jl](https://github.com/JuliaIntervals/IntervalLinearAlgebra.jl) -- Routines to perform numerical linear algebra using interval arithmetic.
- [InvariantSets.jl](https://github.com/ueliwechsler/InvariantSets.jl) -- Compute, approximate and display invariant sets.
- [InvariantSetApproximation.jl](https://github.com/psace-uofa/InvariantSetApproximation.jl) -- Invariant sets of general discrete-time dynamical systems with controls and uncertainties using graph-based algorithms.
- [NeuralVerification.jl](https://github.com/sisl/NeuralVerification.jl) -- Methods to verify deep neural networks.
- [OpticSim.jl](https://github.com/microsoft/OpticSim.jl) -- Ray tracing for procedurally generated systems.
- [Photometry.jl](https://github.com/JuliaAstro/Photometry.jl) -- Utilities for characterizing sources in astronomical images.
- [ReachabilityAnalysis.jl](https://github.com/JuliaReach/ReachabilityAnalysis.jl) -- Methods to compute the sets of states reachable in dynamical systems.
- [Swalbe.jl](https://github.com/Zitzeronion/Swalbe.jl) -- Simple Julia Lattice Boltzmann Solver for Thin Liquid Films and Droplets.
- [TrajectoryGamesBase.jl](https://github.com/lassepe/TrajectoryGamesBase.jl) -- Interface to define trajectory games.

</details>

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


## 📜 How to cite

If you use this package in your work, please cite it using the metadata [here](CITATION.bib) or below.

<details>
<summary>Click to see BibTeX entry. </summary>

```
@article{lazysets21,
  title     = {{LazySets.jl: Scalable Symbolic-Numeric Set Computations}},
  author    = {Forets, Marcelo and Schilling, Christian},
  journal   = {Proceedings of the JuliaCon Conferences},
  year      = {2021},
  publisher = {The Open Journal},
  volume    = {1},
  number    = {1},
  pages     = {11},
  doi       = {10.21105/jcon.00097}
}
```

</details>
