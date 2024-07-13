"""
    σ(d::AbstractVector, P::VPolytope)

Return a support vector of a polytope in vertex representation in a given
direction.

### Input

- `d` -- direction
- `P` -- polytope in vertex representation

### Output

A support vector in the given direction.

### Algorithm

A support vector maximizes the support function.
For a polytope, the support function is always maximized in some vertex.
Hence it is sufficient to check all vertices.
"""
function σ(d::AbstractVector, P::VPolytope)
    return _σ_vertices(d, P.vertices)
end

function _σ_vertices(d, vlist)
    # base cases
    m = length(vlist)
    if m == 0
        error("the support vector of an empty polytope is undefined")
    elseif m == 1
        @inbounds return vlist[1]
    end

    # evaluate support function in every vertex
    N = promote_type(eltype(d), eltype(@inbounds vlist[1]))
    max_ρ = N(-Inf)
    max_idx = 0
    for (i, vi) in enumerate(vlist)
        ρ_i = dot(d, vi)
        if ρ_i > max_ρ
            max_ρ = ρ_i
            max_idx = i
        end
    end
    @inbounds return vlist[max_idx]
end
