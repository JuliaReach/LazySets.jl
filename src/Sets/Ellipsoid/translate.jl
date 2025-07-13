@validate function translate!(E::Ellipsoid, v::AbstractVector)
    c = E.center
    c .+= v
    return E
end
