function ispolyhedral(cp::CartesianProduct)
    ispolyhedraltype(typeof(cp)) && return true

    return ispolyhedral(cp.X) && ispolyhedral(cp.Y)
end
