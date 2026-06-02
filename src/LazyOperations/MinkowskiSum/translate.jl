@validate function translate(ms::MinkowskiSum, x::AbstractVector)
    X = translate(first(ms), x)
    Y = translate(second(ms), x)
    return MinkowskiSum(X, Y)
end
