@validate function translate(P::SimpleSparsePolynomialZonotope, v::AbstractVector)
    return translate!(copy(P), v)
end

function translate!(P::SimpleSparsePolynomialZonotope, v::AbstractVector)
    center(P) .+= v
    return P
end
