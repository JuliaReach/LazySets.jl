"""
    reflect(P::VPolytope)

Concrete reflection of a polytope in vertex representation `P`, resulting in the
reflected set `-P`.

### Input

- `P` -- polytope in vertex representation

### Output

The `VPolytope` representing `-P`.
"""
function reflect(P::VPolytope)
    return VPolytope(-P.vertices)
end
