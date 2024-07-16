function scale!(α::Real, H::Hyperrectangle)
    H.center .*= α
    H.radius .*= abs(α)
    return H
end
