"""
    distance(x::AbstractVector, X::LazySet; [p]::Real=2)
    distance(X::LazySet, x::AbstractVector; [p]::Real=2)

Compute the standard distance (induced by the ``p``-norm) between a point and a
set.

### Input

- `x` -- point/vector
- `X` -- set
- `p` -- (optional; default: `2`) value of the ``p``-norm

### Output

A real number representing the distance between `x` and `X`.

### Notes

The standard distance is zero if the point lies inside the set, and otherwise it
is the ``p``-norm of the shortest line segment between any pair of points.
Formally,

```math
    \\inf_{x' ∈ X} \\{ d(x, x') \\}.
```
"""
@commutative function distance(::AbstractVector, ::LazySet; p::Real=2) end
