function scale(α::Real, H::Hyperrectangle)
    return _scale_copy_inplace(α, H)
end

function scale!(α::Real, H::Hyperrectangle)
    H.center .*= α
    H.radius .*= abs(α)
    return H
end
