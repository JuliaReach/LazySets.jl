"""
    ρ(d::AbstractVector, P::SparsePolynomialZonotope; [enclosure_method]=nothing)

Bound the support function of ``P`` in the direction ``d``.

### Input

- `d`                -- direction
- `P`                -- sparse polynomial zonotope
- `enclosure_method` -- (optional; default: `nothing`) method to use for
                        enclosure; an `AbstractEnclosureAlgorithm` from the
                        [`Rangeenclosures.jl`](https://github.com/JuliaReach/RangeEnclosures.jl)
                        package

### Output

An overapproximation of the support function in the given direction.

### Algorithm

This method implements Proposition 3.1.16 in [1].

[1] Kochdumper, Niklas. *Extensions of polynomial zonotopes and their application to verification of cyber-physical systems.*
    PhD diss., Technische Universität München, 2022.
"""
function ρ(d::AbstractVector, P::SparsePolynomialZonotope;
           enclosure_method=nothing)
    require(@__MODULE__, :RangeEnclosures; fun_name="ρ")
    return _ρ_range_enclosures(d, P, enclosure_method)
end

function load_RangeEnclosures_rho()
    return quote
        function _ρ_range_enclosures(d::AbstractVector, P::SparsePolynomialZonotope,
                                     method::Union{RangeEnclosures.AbstractEnclosureAlgorithm,
                                                   Nothing})
            # default method: BranchAndBoundEnclosure
            isnothing(method) && (method = RangeEnclosures.BranchAndBoundEnclosure())

            c = center(P)
            G = genmat_dep(P)
            GI = genmat_indep(P)
            E = expmat(P)
            n = dim(P)

            res = d' * c + sum(abs.(d' * gi) for gi in eachcol(GI); init=zero(eltype(GI)))

            f(x) = sum(d' * gi * prod(x .^ ei) for (gi, ei) in zip(eachcol(G), eachcol(E)))

            dom = IA.IntervalBox(IA.interval(-1, 1), n)
            res += IA.sup(enclose(f, dom, method))
            return res
        end
    end
end
