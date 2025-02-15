function scale(α::Real, L::LineSegment)
    return _scale_copy_inplace(α, L)
end

function scale!(α::Real, L::LineSegment)
    L.p .*= α
    L.q .*= α
    return L
end
