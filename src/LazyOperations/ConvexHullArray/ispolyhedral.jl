function ispolyhedral(cha::ConvexHullArray)
    ispolyhedraltype(typeof(cha)) && return true

    return all(ispolyhedral, array(cha))
end
