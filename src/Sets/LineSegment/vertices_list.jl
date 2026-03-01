@validate function vertices_list(L::LineSegment)
    return [L.p, L.q]
end
