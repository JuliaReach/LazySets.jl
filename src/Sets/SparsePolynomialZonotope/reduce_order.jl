"""
# Extended help

    reduce_order(P::SparsePolynomialZonotope, r::Real,
                 [method]::AbstractReductionMethod=GIR05())

### Notes

This method implements the algorithm described in [Kochdumper21a; Proposition 3.1.39](@citet).
"""
function reduce_order(P::SparsePolynomialZonotope, r::Real,
                      method::AbstractReductionMethod=GIR05())
    @assert r ≥ 1 "cannot reduce below order 1 (got $r)"

    if order(P) ≤ r
        return P
    end

    n = dim(P)
    h = ngens_dep(P)
    q = ngens_indep(P)

    a = ceil(Int, h + q - n * (r - 1))
    @assert n ≤ a ≤ h + q "unexpected state"  # holds because `1 ≤ r < order(P)`

    G = genmat_dep(P)
    GI = genmat_indep(P)

    # sort all generators by their norm
    Gbar = hcat(G, GI)
    norms = [norm(g) for g in eachcol(Gbar)]
    threshold = sort(norms)[a]

    return _absorb_generators_spz(P, norms, threshold, method)
end

# see ext/LazySets/LazySetsSparsePolynomialZonotopeExt.jl
_absorb_generators_spz(P, norms, threshold, method) = error()
