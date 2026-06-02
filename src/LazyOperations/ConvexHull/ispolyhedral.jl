function ispolyhedral(ch::ConvexHull)
    ispolyhedraltype(typeof(ch)) && return true

    return ispolyhedral(ch.X) && ispolyhedral(ch.Y)
end
