module LazySetsRangeEnclosuresExt

using IntervalArithmetic: interval, sup
using LazySets: AbstractSparsePolynomialZonotope, center, dim, expmat,
                genmat_dep, genmat_indep
using RangeEnclosures: AbstractEnclosureAlgorithm,  # NOTE: this is an internal function
                       BranchAndBoundEnclosure, enclose
import LazySets: _ρ_range_enclosures

function _ρ_range_enclosures(d::AbstractVector, P::AbstractSparsePolynomialZonotope,
                             method::Union{AbstractEnclosureAlgorithm,Nothing})
    # default method: BranchAndBoundEnclosure
    isnothing(method) && (method = BranchAndBoundEnclosure())

    c = center(P)
    G = genmat_dep(P)
    GI = genmat_indep(P)
    E = expmat(P)
    n = dim(P)

    res = d' * c + sum(abs.(d' * gi) for gi in eachcol(GI); init=zero(eltype(GI)))

    f(x) = sum(d' * gi * prod(x .^ ei) for (gi, ei) in zip(eachcol(G), eachcol(E)))

    dom = fill(interval(-1, 1), n)
    res += sup(enclose(f, dom, method))
    return res
end

end  # module
