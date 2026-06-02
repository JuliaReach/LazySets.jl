"""
    in(x::AbstractVector, ia::IntersectionArray)

Check whether a given point is contained in an intersection of a finite number
of sets.

### Input

- `x`  -- point/vector
- `ia` -- intersection of a finite number of sets

### Output

`true` iff ``x ∈ ia``.

### Algorithm

A point ``x`` is in the intersection iff it is in each set.
"""
@validate function in(x::AbstractVector, ia::IntersectionArray)
    for S in ia.array
        if x ∉ S
            return false
        end
    end
    return true
end
