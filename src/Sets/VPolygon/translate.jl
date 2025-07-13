@validate function translate(P::VPolygon, v::AbstractVector)
    return translate!(deepcopy(P), v)
end

function translate!(P::VPolygon, v::AbstractVector)
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    for x in P.vertices
        x .+= v
    end
    return P
end
