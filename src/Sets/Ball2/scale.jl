function scale(α::Real, B::Ball2)
    return Ball2(B.center .* α, B.radius * abs(α))
end
