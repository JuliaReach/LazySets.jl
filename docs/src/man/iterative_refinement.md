# Iterative Refinement

This section of the manual describes an approximation method for an arbitrary
two-dimensional convex set ``S`` and a given error bound ``ɛ`` using support
vectors.

```@meta
CurrentModule = LazySets.Approximations
```

The basic idea is to add new supporting directions whenever the approximation
error is still bigger than ``ɛ``.

The approximation is represented by a list of local refinements.
Each refinement describes a set with one angle and is wrapped in the following
type.

The approximation is initialized with box directions, i.e., we have four
refinement instances, one for each angle.
Then we just iterate through all refinement instances and check if the error is
bigger than the threshold individually.
If so, we refine the instance by splitting into two more precise refinement
instances and apply the checks recursively.
