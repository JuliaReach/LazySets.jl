"""
    ρ(d::AbstractVector, P::VPolytope)

Evaluate the support function of a polytope in vertex representation in a given
direction.

### Input

- `d` -- direction
- `P` -- polytope in vertex representation

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector, P::VPolytope)
    return _ρ_vertices(d, P.vertices)
end

function _ρ_vertices(d, vlist)
    if isempty(vlist)
        error("the support function of an empty polytope is undefined")
    end
    # evaluate support function in every vertex
    return maximum(v -> dot(d, v), vlist)
end
