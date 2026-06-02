@validate function translate(ia::IntersectionArray, x::AbstractVector)
    return IntersectionArray([translate(X, x) for X in array(ia)])
end
