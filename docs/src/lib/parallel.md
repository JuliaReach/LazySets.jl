# Parallel

This section of the manual describes the parallel implementation of some
algorithms in the `LazySets` library. These algorithms are implemented in the
`LazySets.Parallel` module.

```@contents
Pages = ["parallel.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Parallel
DocTestSetup = quote
    using LazySets, LazySets.Parallel
end
```

```@docs
Parallel
```

## Box approximations

```@docs
LazySets.Parallel.interval_hull
LazySets.Parallel.box_approximation
LazySets.Parallel.symmetric_interval_hull
LazySets.Parallel.box_approximation_symmetric
LazySets.Parallel.box_approximation_helper_parallel
LazySets.Parallel.process_chunk!
```

## Distributed functions

```@docs
LazySets.Parallel.assign_chunk!
LazySets.Parallel.distribute_task!
LazySets.Parallel.prange
```
