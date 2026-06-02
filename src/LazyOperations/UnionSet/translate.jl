@validate function translate(cup::UnionSet, x::AbstractVector)
    X = translate(first(cup), x)
    Y = translate(second(cup), x)
    return UnionSet(X, Y)
end
