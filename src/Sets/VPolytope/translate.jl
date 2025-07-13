@validate function translate(P::VPolytope, v::AbstractVector)
    return translate!(deepcopy(P), v)
end

function translate!(P::VPolytope, v::AbstractVector)
    if isempty(P.vertices)
        return P
    end

    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    for x in P.vertices
        x .+= v
    end
    return P
end
