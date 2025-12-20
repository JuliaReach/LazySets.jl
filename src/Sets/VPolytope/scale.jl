function scale(α::Real, P::VPolytope)
    return VPolytope([α * v for v in P.vertices])
end

function scale!(α::Real, P::VPolytope)
    @inbounds for v in P.vertices
        v .*= α
    end
    return P
end
