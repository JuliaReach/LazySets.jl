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
to launch as many workers as the number of local CPU threads.


```@contents
Pages = ["parallel_approximations.md"]
Depth = 3
```

```@meta
DocTestSetup = quote
    using LazySets, LazySets.Parallel
end
```

## Parallel interval hulls 

Consider the symmetric interval hull approximation
of a nested lazy set in `n` dimensions, where `n` is a positive integer.
The calculation, a modified version of that given in the Introduction,
and that is used for example in the discretization of set-based ODEs, is defined
below.

```julia
using LazySets
using SparseArrays

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
benchmarking. We will explore the increase in the computational cost by increasing
the dimension `n`, and compare the sequential algorithm with the parallel algorithm.

Given the lazy set `Y(n)`, we want to calculate the symmetric interval hull, which
corresponds to finding the smallest `n`-dimensional hyperrectangle that contains
the set `Y(n)` and is symmetric with respect to the origin. Notice that this operation
is inherently parallel, since one can evaluate the support function of `Y` independently
on each dimension from `1` to `n`.

The sequential algorithm returns the following execution times. We use
the `@btime` macro from the `BenchmarkTools` package to have a more accurate
timing than `@time`; the `$n` argument is used for interpolation of the arguments
(if you are not behchmarking, pass `n` to `symmetric_interval_hull`, as usual).

```julia
using BenchmarkTools

for n in [50, 100, 500, 1000]
    @btime res = Approximations.symmetric_interval_hull(Y($n));
end

  60.690 ms (11653 allocations: 25.89 MiB)
  140.726 ms (23317 allocations: 54.16 MiB)
  2.211 s (116529 allocations: 381.27 MiB)
  11.527 s (234505 allocations: 1.01 GiB)
```

For the parallel benchmark, we start Julia with 4 processes with the command
`$ julia -p 4` and call `LazySets.Parallel.symmetric_interval_hull(Y(n))`. 

```julia
import LazySets.Parallel

for n in [50, 100, 500, 1000]
    @time LazySets.Parallel.symmetric_interval_hull(Y(n));
end

  0.008020 seconds (2.89 k allocations: 170.063 KiB)
  0.016781 seconds (3.63 k allocations: 273.453 KiB)
  0.361447 seconds (11.24 k allocations: 2.505 MiB)
  2.486827 seconds (22.24 k allocations: 8.762 MiB)
```

In the following table we summarize the speedup.

|n|Sequential (s)| Parallel `p=4` (s) | Speedup|
|---|----|----|----|
|50| 0.061  | 0.008 | 7.62|
|100| 0.14 | 0.016 | 8.75 |
|500| 2.2  | 0.36 | 6.11|
|1000| 11.5 | 2.48 | 4.63|

We haven't used `@btime` macro for the parallel benchmark, because it gives
`ERROR: SystemError: shm_open() ... Too many open files`.

!!! note
    The results in this section were obtained with a standard MacBook Pro laptop
    with the following specifications:
    ```julia
    Julia Version 0.7.0
    Commit a4cb80f3ed (2018-08-08 06:46 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin14.5.0)
      CPU: Intel(R) Core(TM) i7-4770HQ CPU @ 2.20GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-6.0.0 (ORCJIT, haswell)
    ```
