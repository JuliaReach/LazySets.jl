@validate function translate(P::SparsePolynomialZonotope, v::AbstractVector)
    return translate!(copy(P), v)
end

@validate function translate!(P::SparsePolynomialZonotope, v::AbstractVector)
    center(P) .+= v
    return P
end
