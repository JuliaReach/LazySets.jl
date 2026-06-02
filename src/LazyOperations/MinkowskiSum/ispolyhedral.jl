function ispolyhedral(ms::MinkowskiSum)
    ispolyhedraltype(typeof(ms)) && return true

    return ispolyhedral(ms.X) && ispolyhedral(ms.Y)
end
