# Parallel Approximations

A subset of the approximation algorithms are implemented in parallel in the 
`LazySets.Parallel` module. In order to use parallel versions of the algorithms,
you can write:

```@example
using LazySets
import LazySets.Parallel

# call a method implemented in parallel, for example:
S = Ball2(ones(100), 1.0)
Parallel.box_approximation(S)
```

Note that after importing or using `LazySets.Parallel`, the version of the function
used must be fully qualified, eg. `LazySets.Approximations.box_approximation` for the
sequential version or `LazySets.Parallel.box_approximation` for the parallel version.

The parallelization strategy that is available uses processes. To set the number
of processes `N`, use the flag `-p N` at julia startup. For example, do

```julia
$ julia -p 4
```
to launch `4` additional local worker julia processes. Use the keyword `auto`,
as in
```julia
$ julia -p auto
```
to launch as many workers as the number of local CPU cores.


```@contents
Pages = ["parallel_approximations.md"]
Depth = 3
```

## Parallel interval hulls 

As an illustration of the symmetric interval hull approximation of a nested
lazy set computed in parallel, consider the following calculation.
It arises in the discretization of set-based ODEs, and is defined below for an artificial
example of a tridiagonal matrix of order `n`, where `n` is a positive integer.

```julia
using LazySets, Expokit
using SparseArrays, LinearAlgebra

# define an nxn tridiagonal matrix
A(n) = sparse(diagm(0 => fill(0.05, n), -1 => fill(-1, n-1), 1 => fill(-1, n-1)))

# step size and initial set
δ = 0.1
X0(n) = Ball2(ones(n), 0.1)

# input coefficients matrix (nx2 matrix with coefficients from -1 to 1)
b(n) = vcat(range(-1, stop=1, length=n))
B(n) = [b(n) b(n)] 
U = BallInf(zeros(2), 1.2)

# lazy matrix exponential
eAδ(n) = SparseMatrixExp(A(n) * δ)

# set that we want to overapproximate with an interval hull
Y(n) = ConvexHull(eAδ(n) * X0(n) ⊕ (δ * B(n) * U), X0(n))
```

The set `Y(n)` is parametric in the system's dimension `n`, to facilitate
benchmarking. We will explore the computational cost as the dimension `n` increases,
and compare the sequential algorithm with the parallel algorithm.

Given the lazy set `Y(n)`, we want to calculate the symmetric interval hull, which
corresponds to finding the smallest `n`-dimensional hyperrectangle that contains
the set `Y(n)` and is symmetric with respect to the origin. Notice that this operation
is inherently parallel, since one can evaluate the support function of `Y` independently
in each dimension from `1` to `n`.

The sequential algorithm returns the following execution times. We use
the `@btime` macro from the `BenchmarkTools` package to have a more accurate
timing than `@time`; the `$n` argument is used for interpolation of the arguments
(if you are not behchmarking, pass `n` to `symmetric_interval_hull`, as usual).

```julia
using BenchmarkTools

for n in [50, 100, 500, 1000]
    @btime res = Approximations.symmetric_interval_hull(Y($n));
end

  59.103 ms (11554 allocations: 25.89 MiB)
  129.453 ms (23118 allocations: 54.16 MiB)
  1.943 s (115530 allocations: 381.26 MiB)
  10.017 s (232506 allocations: 1.01 GiB)
```

For the parallel benchmark, we start Julia with 4 processes with the command
`$ julia -p 4` and call `LazySets.Parallel.symmetric_interval_hull(Y(n))`. 

```julia
import LazySets.Parallel

for n in [50, 100, 500, 1000]
    @btime LazySets.Parallel.symmetric_interval_hull($Y($n));
end

  6.846 ms (2550 allocations: 160.59 KiB)
  13.544 ms (3528 allocations: 271.94 KiB)
  387.556 ms (11155 allocations: 2.51 MiB)
  2.638 s (22156 allocations: 8.77 MiB)
```

In the following table we summarize the speedup.

|n|Sequential (s)| Parallel `p=4` (s) | Speedup|
|---|----|----|----|
|50| 0.059  | 0.007 | 8.42|
|100| 0.129 | 0.013 | 9.92 |
|500| 1.94  | 0.387 | 4.96|
|1000| 10.0 | 2.64 | 3.79|

The results in this section were obtained with a standard MacBook Pro laptop
with the following specifications:

```julia
julia> versioninfo()
Julia Version 1.0.2
Commit d789231e99 (2018-11-08 20:11 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin14.5.0)
  CPU: Intel(R) Core(TM) i7-4770HQ CPU @ 2.20GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.0 (ORCJIT, haswell)
```
