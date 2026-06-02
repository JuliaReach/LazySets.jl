@validate function translate(cap::Intersection, x::AbstractVector)
    X = translate(first(cap), x)
    Y = translate(second(cap), x)
    return Intersection(X, Y; cache=cap.cache)
end
