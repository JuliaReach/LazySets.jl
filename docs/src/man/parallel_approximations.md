# Parallel Approximations

A subset of the approximation algorithms are implemented in parallel in the 
`LazySets.Parallel` module. In order to use parallel versions of the algorithms,
you can write:

```@example
using LazySets
import LazySets.Parallel

S = Ball2(ones(100), 1.0)
Parallel.box_approximation(S)
```

Note that after importing or using `LazySets.Parallel`, the version of the function
used must be fully qualified, eg. `LazySets.Approximations.box_approximation` for the
sequential function or `LazySets.Parallel.box_approximation` for the parallel version).

---

```@contents
Pages = ["parallel_approximations.md"]
Depth = 3
```

```@meta
DocTestSetup = quote
    using LazySets, LazySets.Parallel
end
```