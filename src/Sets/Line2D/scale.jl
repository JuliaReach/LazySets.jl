function scale(α::Real, L::Line2D)
    return Line2D(copy(L.a), α * L.b)
end

function scale!(α::Real, L::Line2D)
    L.a ./= α
end
