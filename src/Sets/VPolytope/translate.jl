@validate function translate(P::VPolytope, v::AbstractVector)
    return translate!(deepcopy(P), v)
end

@validate function translate!(P::VPolytope, v::AbstractVector)
    if isempty(P.vertices)
        return P
    end

    for x in P.vertices
        x .+= v
    end
    return P
end
