function scale!(α::Real, L::LineSegment)
    L.p .*= α
    L.q .*= α
    return L
end
