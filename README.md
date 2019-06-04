# LazySets.jl

[![Build Status](https://travis-ci.org/JuliaReach/LazySets.jl.svg?branch=master)](https://travis-ci.org/JuliaReach/LazySets.jl)
[![Docs latest](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliareach.github.io/LazySets.jl/latest/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/JuliaReach/LazySets.jl/blob/master/LICENSE.md)
[![DOI](https://zenodo.org/badge/105701832.svg)](https://zenodo.org/badge/latestdoi/105701832)
[![Code coverage](http://codecov.io/github/JuliaReach/LazySets.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaReach/LazySets.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

`LazySets` is a [Julia](http://julialang.org) package for calculus with convex sets.

## Resources

- [Manual](http://juliareach.github.io/LazySets.jl/latest/)
- [Contributing](https://juliareach.github.io/LazySets.jl/latest/about.html#Contributing-1)
- [Release notes of tagged versions](https://github.com/JuliaReach/LazySets.jl/releases)
- [Release notes of the development version](https://github.com/JuliaReach/LazySets.jl/wiki/Release-log-tracker)
- [Publications](https://juliareach.github.io/Reachability.jl/latest/publications.html)

## Installing

This package requires Julia v1.0 or later.
Refer to the [official documentation](https://julialang.org/downloads) on how to
install and run Julia on your system.

Depending on your needs, choose an appropriate command from the following list
and enter it in Julia's REPL.
To activate the `pkg` mode, type `]` (and to leave it, type `exit()`).

#### [Install the latest release version](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Adding-registered-packages-1)

```julia
pkg> add LazySets
```

#### Install the latest development version

```julia
pkg> add LazySets#master
```

#### [Clone the package for development](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Developing-packages-1)

```julia
pkg> dev LazySets
```
