function scale(α::Real, P::VPolygon)
    return VPolygon([α * v for v in P.vertices])
end

function scale!(α::Real, P::VPolygon)
    @inbounds for v in P.vertices
        v .*= α
    end
    return P
end
