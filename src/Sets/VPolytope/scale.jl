function scale(α::Real, P::VPolytope)
    return _scale_copy_inplace(α, P)
end

function scale!(α::Real, P::VPolytope)
    P.vertices .*= α
    return P
end
