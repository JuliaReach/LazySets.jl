@validate function translate(cp::CartesianProduct, x::AbstractVector)
    X = first(cp)
    n1 = dim(X)
    X = translate(X, @view(x[1:n1]))
    Y = translate(second(cp), @view(x[(n1 + 1):end]))
    return CartesianProduct(X, Y)
end
