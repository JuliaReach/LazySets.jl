function reduce_order(A::MatrixZonotope, r::Real,
                      method::AbstractReductionMethod=GIR05())
    @assert r ≥ 1 "cannot reduce below order 1 (got $r)"

    if order(A) <= r
        return A
    end

    Z = convert(Zonotope, A)
    Zred = reduce_order(Z, r, method)

    # reshape to matrix zonotope
    gens = [Matrix(reshape(col, size(A))) for col in eachcol(genmat(Zred))]
    return MatrixZonotope(center(A), gens)
end