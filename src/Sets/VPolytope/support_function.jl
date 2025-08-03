@validate function ρ(d::AbstractVector, P::VPolytope)
    return _ρ_vertices(d, P.vertices)
end

function _ρ_vertices(d, vlist)
    # @validate ensures `!isempty(vlist)`

    # evaluate support function in every vertex
    return maximum(v -> dot(d, v), vlist)
end
