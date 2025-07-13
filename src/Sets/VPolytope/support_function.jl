@validate function ρ(d::AbstractVector, P::VPolytope)
    return _ρ_vertices(d, P.vertices)
end

function _ρ_vertices(d, vlist)
    if isempty(vlist)
        error("the support function of an empty polytope is undefined")
    end
    # evaluate support function in every vertex
    return maximum(v -> dot(d, v), vlist)
end
