function scale(α::Real, P::Polygon)
    return Polygon([α * v for v in P.vertices])
end

function scale!(α::Real, P::Polygon)
    @inbounds for v in P.vertices
        v .*= α
    end
    return P
end
