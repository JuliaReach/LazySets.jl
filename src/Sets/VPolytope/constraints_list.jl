"""
    constraints_list(P::VPolytope)

Return a list of constraints defining a polytope in vertex representation.

### Input

- `P` -- polytope in vertex representation

### Output

A list of constraints of the polytope.

### Algorithm

We use `tohrep` to compute the constraint representation of `P`.
"""
function constraints_list(P::VPolytope)
    n = dim(P)
    if n == 1
        return constraints_list(convert(Interval, P))
    elseif n == 2
        return constraints_list(convert(VPolygon, P))
    else
        return constraints_list(tohrep(P))
    end
end
