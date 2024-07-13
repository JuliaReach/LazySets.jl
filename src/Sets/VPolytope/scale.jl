function scale!(α::Real, P::VPolytope)
    P.vertices .*= α
    return P
end
