function ispolyhedral(cap::Intersection)
    ispolyhedraltype(typeof(cap)) && return true

    return ispolyhedral(cap.X) && ispolyhedral(cap.Y)
end
