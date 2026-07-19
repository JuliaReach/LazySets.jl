using IntervalArithmetic: isequal_interval, mid, radius
using LazySets: center, dim, expmat, genmat_dep, genmat_indep, ngens_dep,
                ngens_indep, nparams, polynomial_order, sym_box, sym_itv,
                zero_box, zero_itv
using LazySets.SparsePolynomialZonotopeModule: SparsePolynomialZonotope
using LinearAlgebra: diagm
using ReachabilityBase.Arrays: SingleEntryVector, remove_zero_columns
using TaylorModels: TaylorModelN, constant_term, polynomial, remainder
using TaylorModels.TaylorSeries: HomogeneousPolynomial, TaylorN, set_variables,
                                 coeff_table, in_base, pos_table  # NOTE: these are internal functions
import Base: convert

# check that a vector of Taylor models has the [-1, 1] domain
function _has_normalized_domain(vTM::Vector)
    return all(TMi -> all(isequal_interval(sym_itv(_eltype_TM(TMi))), domain(TMi)), vTM)
end

# implements Proposition 3.1.12 in thesis
function convert(::Type{SparsePolynomialZonotope}, vTM::Vector{<:TaylorModelN{r,N}}) where {r,N}
    @assert _has_normalized_domain(vTM) "normalized domain (-1, 1) required"

    # upper bound on the number of terms/columns (ignores duplicates)
    # - 1 per iteration because we treat the constant terms separately
    num_coeffs = 0
    for TMi in vTM
        sum(TMi -> length(polynomial(TMi).coeffs) - 1, vTM)
    end

    n = length(vTM)
    c = Vector{N}(undef, n)
    Gs = Vector{Vector{N}}()
    GI_diag = Vector{N}(undef, n)
    Es = Vector{Vector{Int}}()

    total_columns = 0
    @inbounds for (i, TMi) in enumerate(vTM)
        pol = polynomial(TMi)
        rem = remainder(TMi)
        c[i] = mid(rem) + constant_term(pol)
        GI_diag[i] = radius(rem)
        for (order, term_order) in enumerate(pol.coeffs)
            # we skip the first (= constant) term
            if order == 1 || iszero(term_order)
                continue
            end
            for (k, coeff_k) in enumerate(term_order.coeffs)
                if iszero(coeff_k)
                    continue
                end
                Ej = coeff_table[order][k]
                j = findfirst(e -> e == Ej, Es)
                if isnothing(j)
                    total_columns += 1
                    j = total_columns
                    push!(Es, Ej)
                    push!(Gs, zeros(N, n))
                end
                Gs[j][i] += coeff_k
            end
        end
    end
    G = stack(Gs)
    GI = remove_zero_columns(diagm(GI_diag))
    E = stack(Es)
    return SparsePolynomialZonotope(c, G, GI, E)
end

# implements Proposition 3.1.13 in thesis
function convert(::Type{Vector{<:TaylorModelN}}, P::SparsePolynomialZonotope{N}) where {N}
    n = dim(P)
    p = nparams(P)  # number of parameters
    q = ngens_indep(P)  # independent/linear parameters
    h = ngens_dep(P)  # dependent/nonlinear parameters
    r = p + q
    c = center(P)
    G = genmat_dep(P)
    GI = genmat_indep(P)
    E = expmat(P)
    poly_order = polynomial_order(P)
    z = zeros(Int, q)
    # we need to rewrite the global variables
    set_variables("x"; order=poly_order, numvars=r)
    # the following vectors are shared for each polynomial
    rem = zero_itv(N)
    dom = sym_box(r, N)
    x0 = zero_box(r, N)

    vTM = Vector{TaylorModelN{r,N,N}}(undef, n)
    @inbounds for i in 1:n
        coeffs = zeros(HomogeneousPolynomial{N}, poly_order)

        # constant term
        coeffs[1] = c[i]

        # linear terms
        for j in 1:q
            v = zeros(N, r)
            v[p + j] = GI[i, j]
            coeffs[2] += HomogeneousPolynomial(v, 1)
        end

        # nonlinear terms
        # TODO can be faster if loop is restructured so that index j is
        # created for each dimension
        for j in 1:h
            Ej = E[:, j]
            ord = sum(Ej)
            G[i, j]
            idx = pos_table[ord + 1][in_base(poly_order, Ej)]
            l = length(coeff_table[ord + 1])
            v = Vector(SingleEntryVector(idx, l, G[i, j]))
            pj = HomogeneousPolynomial(v, ord)
            coeffs[ord + 1] += pj
        end

        pol = TaylorN(coeffs, poly_order)
        vTM[i] = TaylorModelN(pol, rem, x0, dom)
    end
    return vTM
end
