function center(L::LineSegment)
    return L.p + (L.q - L.p) / 2
end
