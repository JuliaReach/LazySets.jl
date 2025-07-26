"""
	linear_map(MZ::MatrixZonotope, P::SparsePolynomialZonotope)

Compute the linear map of a sparse polynomial zonotope through a matrix zonotope,
following Proposition 1 of [HuangLBS2025](@citet).

### Input

- `lm` -- a linear map of a sparse polynomial zonotope through a matrix zonotope

### Output

A sparse polynomial zonotope representing the linear map ``MZ ⋅ P```.

"""
function linear_map(MZ::MatrixZonotope, P::SparsePolynomialZonotope)
    if ngens_indep(P) > 0
        error("An exact expression for the linear map is only available for " *
              "`SparsePolynomialZonotope`s with no independent generators. " *
              "Try using `overapproximate` instead.")
    end

    @assert size(MZ, 2) == dim(P) "a linear map of size $(size(M)) cannot " *
                                  "be applied to a set of dimension $(dim(X))"

    T = promote_type(eltype(MZ), eltype(P))

    n = dim(P)
    w = ngens(MZ)
    h = ngens_dep(P)

    c = center(MZ) * center(P)

    # compute matrix of dependent generators
    G = Matrix{T}(undef, n, h + w + h * w)
    G[:, 1:h] = center(MZ) * genmat_dep(P)
    @inbounds for (i, A) in enumerate(generators(MZ))
        G[:, h + i] = A * center(P)
        G[:, (h + w + (i - 1) * h + 1):(h + w + i * h)] = A * genmat_dep(P)
    end

    Gi = Matrix{T}(undef, n, 0)

    # compute exponent
    Imat = Matrix{Int}(I, w, w)
    Ê₁, Ê₂, idx = merge_id(indexvector(P), indexvector(MZ), expmat(P), Imat)
    pₖ = size(Ê₁, 1)
    E = Matrix{Int}(undef, pₖ, h + w + h * w)
    E[:, 1:h] = Ê₁
    E[:, (h + 1):(h + w)] = Ê₂
    ones_h = ones(Int, 1, h)
    @inbounds for l in 1:w
        cstart = (h + w) + (l - 1) * h + 1
        cend = (h + w) + l * h
        col = @view Ê₂[:, l]
        E[:, cstart:cend] = col * ones_h .+ Ê₁
    end

    return SparsePolynomialZonotope(c, G, Gi, E, idx)
end

"""
	linear_map(MZP::MatrixZonotopeProduct, P::SparsePolynomialZonotope)

Compute the linear map of a sparse polynomial zonotope through a matrix zonotope product
by recursively applying the overapproximation rule from the inside out.

### Input

- `MZ` -- a matrix zonotope product
- `P` -- a sparse polynomial zonotope 

### Output

A sparse polynomial zonotope representing the linear map ``MZ ⋅ P```.

"""
function linear_map(MZP::MatrixZonotopeProduct, P::SparsePolynomialZonotope)
    MZs = factors(MZP)
    reduced = foldr((A, acc) -> linear_map(A, acc), MZs; init=P)
    return reduced
end
