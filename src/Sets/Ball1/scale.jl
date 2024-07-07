function scale(α::Real, B::Ball1)
    return Ball1(B.center .* α, B.radius * abs(α))
end
