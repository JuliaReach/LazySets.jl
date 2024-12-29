"""
    distance(x::AbstractVector, X::LazySet; [p]::Real=2)
    distance(X::LazySet, x::AbstractVector; [p]::Real=2)

Compute the standard distance (induced by the ``p``-norm) between a point and a set.

### Input

- `x` -- point/vector
- `X` -- set
- `p` -- (optional; default: `2`) value of the ``p``-norm

### Output

A real number representing the distance between `x` and `X`.

### Notes

The standard distance is zero if the point lies inside the set, and infinite if the set is empty.
Otherwise, it is the ``p``-norm of the shortest line segment between the point and any other point
in the set. Formally,

```math
    \\inf_{y ∈ X} \\{ d(x, y) \\}.
```
"""
function distance(::AbstractVector, ::LazySet; p::Real=2) end
