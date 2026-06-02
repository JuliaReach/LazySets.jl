@validate function translate(cha::ConvexHullArray, x::AbstractVector)
    return ConvexHullArray([translate(X, x) for X in array(cha)])
end
