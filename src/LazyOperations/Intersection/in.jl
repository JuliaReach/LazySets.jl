"""
    in(x::AbstractVector, cap::Intersection)

Check whether a given point is contained in the intersection of two sets.

### Input

- `x`   -- point/vector
- `cap` -- intersection of two sets

### Output

`true` iff ``x ∈ cap``.

### Algorithm

A point ``x`` is in the intersection iff it is in each set.
"""
@validate function in(x::AbstractVector, cap::Intersection)
    return (x ∈ cap.X) && (x ∈ cap.Y)
end
