function scale(α::Real, H::Hyperrectangle)
    return Hyperrectangle(α * H.center, abs(α) * H.radius)
end

function scale!(α::Real, H::Hyperrectangle)
    H.center .*= α
    H.radius .*= abs(α)
    return H
end
