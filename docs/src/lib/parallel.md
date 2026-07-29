# LazySetsDistributedExt

This section of the manual describes the parallel implementation of some
algorithms in the `LazySets` library. These algorithms are implemented in the
`LazySetsDistributedExt` module.

```@contents
Pages = ["parallel.md"]
Depth = 3
```

```@meta
CurrentModule = LazySetsDistributedExt
```

## Box approximations

```@docs
LazySetsDistributedExt.box_approximation
LazySetsDistributedExt.box_approximation_symmetric
LazySetsDistributedExt.box_approximation_helper_parallel
LazySetsDistributedExt.process_chunk!
```

## Distributed functions

```@docs
LazySetsDistributedExt.assign_chunk!
LazySetsDistributedExt.distribute_task!
LazySetsDistributedExt._prange
```
