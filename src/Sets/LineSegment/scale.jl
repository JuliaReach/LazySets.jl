function scale(α::Real, L::LineSegment)
    return LineSegment(α * L.p, α * L.q)
end

function scale!(α::Real, L::LineSegment)
    L.p .*= α
    L.q .*= α
    return L
end
