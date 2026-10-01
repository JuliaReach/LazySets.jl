@validate function translate(P::Polygon, v::AbstractVector)
    return translate!(deepcopy(P), v)
end

@validate function translate!(P::Polygon, v::AbstractVector)
    for x in P.vertices
        x .+= v
    end
    return P
end
