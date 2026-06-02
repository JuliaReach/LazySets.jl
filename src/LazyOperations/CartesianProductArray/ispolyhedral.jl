function ispolyhedral(cpa::CartesianProductArray)
    ispolyhedraltype(typeof(cpa)) && return true

    return all(ispolyhedral, array(cpa))
end
