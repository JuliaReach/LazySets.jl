# Iterative Refinement

This section of the manual describes an approximation method for an arbitrary
two-dimensional convex set ``S`` and a given error bound ``ɛ`` using support
vectors.

```@contents
Pages = ["iterative_refinement.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Approximations
```

The basic idea is to add new supporting directions whenever the approximation
error is still bigger than ``ɛ``.

## Local approximations

The approximation is represented by a list of local approximations or refinements.
Each refinement describes a set with one angle.
More precisely, a *local approximation* is a triple ``(p_1, p_2, q)``, where:

- ``p_1`` and ``p_2`` belong to ``S``
- the segments ``(p_1 q)`` and ``(p_2 q)`` belong to support lines of ``S``

Recall that, since ``S`` is assumed to be convex, the segment ``(p_1 p_2)`` is
inside ``S``. The type `Approximation2D` implements a local approximation in 2D.

## Initialization

The approximation is initialized with box directions, i.e., we have four
refinement instances, one for each angle.

Let `ndir` denote the direction normal to the inner approximation.

## Iteration

Then we just iterate through all refinement instances and check if the error is
bigger than the threshold individually.
If so, we refine the instance by splitting into two more precise refinement
instances and apply the checks recursively.

## Stopping criteria

## Examples

Consider in two dimensions a centered unit ball in the 2-norm. Let's consider a
polygonal approximation with error bound ``ɛ = .01``.

```@example ir
using LazySets, LazySets.Approximations, Plots

b = Ball2(zeros(2), 1.)
ɛ = .5
p = overapproximate(b, ɛ)

plot(b, 1e-3, aspectratio=1)  # use 1e-3 for precise definition of borders
plot!(p, alpha=0.4)
```

Let's analyze in some detail how `p` was actually computed. The "box" directions
are used as the initialization.

```@example ir
import LazySets.Approximations:Approximation2D, refine

const DIR_EAST, DIR_NORTH, DIR_WEST, DIR_SOUTH = [1., 0.], [0., 1.], [-1., 0.], [0., -1.]

```
