"""
# Extended help

    σ(d::AbstractVector, P::VPolytope)

### Algorithm

A support vector maximizes the support function.
For a polytope, the support function is always maximized in some vertex.
Hence it is sufficient to check all vertices.
"""
@validate function σ(d::AbstractVector, P::VPolytope)
    return _σ_vertices(d, P.vertices)
end

function _σ_vertices(d, vlist)
    m = length(vlist)
    # @validate ensures `m > 0`

    if m == 1
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
