@validate function translate(P::SparsePolynomialZonotope, v::AbstractVector)
    return translate!(copy(P), v)
end

function translate!(P::SparsePolynomialZonotope, v::AbstractVector)
    center(P) .+= v
    return P
end
