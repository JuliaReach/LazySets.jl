"""
    isbounded(ia::IntersectionArray)

Check whether an intersection of a finite number of sets is bounded.

### Input

- `ia` -- intersection of a finite number of sets

### Output

`true` iff the intersection is bounded.

### Algorithm

We first check if any of the wrapped sets is bounded.
Otherwise we check boundedness via
[`LazySets._isbounded_unit_dimensions`](@ref).
"""
function isbounded(ia::IntersectionArray)
    if any(isbounded, ia.array)
        return true
    end
    return _isbounded_unit_dimensions(ia)
end
