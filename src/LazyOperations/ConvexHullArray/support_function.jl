"""
    ρ(d::AbstractVector, cha::ConvexHullArray)

Evaluate the support function of a convex hull of a finite number of sets in a
given direction.

### Input

- `d`   -- direction
- `cha` -- convex hull of a finite number of sets

### Output

The evaluation of the support function of the convex hull of a finite number of
sets in the given direction.

### Algorithm

This algorithm calculates the maximum over all ``ρ(d, X_i)``, where the
``X_1, …, X_k`` are the sets in the array of `cha`.
"""
@validate function ρ(d::AbstractVector, cha::ConvexHullArray)
    return maximum(ρ(d, Xi) for Xi in cha)
end
