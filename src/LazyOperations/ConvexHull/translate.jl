@validate function translate(ch::ConvexHull, x::AbstractVector)
    X = translate(first(ch), x)
    Y = translate(second(ch), x)
    return ConvexHull(X, Y)
end
