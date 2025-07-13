@validate function ρ(d::AbstractVector, L::LineSegment)
    return max(dot(L.p, d), dot(L.q, d))
end
