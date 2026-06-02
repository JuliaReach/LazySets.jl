function ispolyhedral(ia::IntersectionArray)
    ispolyhedraltype(typeof(ia)) && return true

    return all(ispolyhedral, array(ia))
end
