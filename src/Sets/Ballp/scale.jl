function scale(α::Real, B::Ballp)
    return Ballp(B.p, B.center .* α, B.radius * abs(α))
end
