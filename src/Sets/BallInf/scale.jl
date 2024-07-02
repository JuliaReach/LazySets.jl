function scale(α::Real, B::BallInf)
    return BallInf(B.center .* α, B.radius * abs(α))
end
