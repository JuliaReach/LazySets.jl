using IntervalArithmetic: isequal_interval, mid, radius
using LazySets: center, dim, expmat, genmat_dep, genmat_indep, ngens_dep,
                ngens_indep, nparams, polynomial_order, sym_box, sym_itv,
                zero_box, zero_itv
using LazySets.SparsePolynomialZonotopeModule: SparsePolynomialZonotope
using LinearAlgebra: diagm
using ReachabilityBase.Arrays: SingleEntryVector, remove_zero_columns
using TaylorModels: TaylorModelN, polynomial, remainder
using TaylorSeries: HomogeneousPolynomial, TaylorN, constant_term, variables!,
                    default_space, in_base  # NOTE: these are internal functions
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
                Ej = pol.space.coeff_table[order][k]
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
    @show P
    @show n = dim(P)
    p = nparams(P)  # number of parameters
    q = ngens_indep(P)  # independent/linear parameters
    @show h = ngens_dep(P)  # dependent/nonlinear parameters
    @show r = p + q
    c = center(P)
    G = genmat_dep(P)
    GI = genmat_indep(P)
    E = expmat(P)
    @show poly_order = polynomial_order(P)
    # we need to rewrite the global variables
    @show space = default_space[]
    if space.num_vars != r || space.order != poly_order
        variables!("x"; numvars=r, order=poly_order, nowarn=true)
        @show 1
        @show space
        @show space.pos_table
    end
    # the following vectors are shared for each polynomial
    rem = zero_itv(N)
    dom = sym_box(r, N)
    x0 = zero_box(r, N)

    vTM = Vector{TaylorModelN{r,N,N}}(undef, n)
    for i in 1:n
        # `coeffs` stores polynomial terms; index i stores terms of order i+1
        # (e.g., index 1 = constant term, index 2 = linear terms, etc.)
        coeffs = zeros(HomogeneousPolynomial{N}, poly_order + 1)

        # constant term from c
        coeffs[1] = c[i]

        # linear terms from GI
        ord = 1
        c_idx = 2
        for j in 1:q
            v = zeros(N, r)
            v[p + j] = GI[i, j]
            pj = HomogeneousPolynomial(space, v, ord)
            if iszero(coeffs[c_idx])
                coeffs[c_idx] = pj
            else
                coeffs[c_idx] += pj
            end
        end
	@show coeffs[2].space

        # linear and nonlinear terms from G & E
        # TODO can be faster if loop is restructured so that index j is
        # created for each dimension
        z_extend = zeros(Int, q)
        for j in 1:h
            @show Ej = E[:, j]
            @show ord = sum(Ej)
            @show c_idx = ord + 1
#             space = coeffs[c_idx].space
            @show space
            @show length(space.pos_table)
            pos_table_entry = space.pos_table[c_idx]
            @show Ej_ext = vcat(Ej, z_extend)
            @show ib = in_base(r, Ej_ext)  # TODO `r` should be `poly_order`, but it does not work
            idx = pos_table_entry[ib]
            l = length(space.coeff_table[c_idx])
            v = Vector(SingleEntryVector(idx, l, G[i, j]))
            @show pj = HomogeneousPolynomial(space, v, ord)
            @show coeffs[c_idx]
            if iszero(coeffs[c_idx])
                coeffs[c_idx] = pj
            else
                coeffs[c_idx] += pj
            end
        end

        pol = TaylorN(coeffs, poly_order)
        vTM[i] = TaylorModelN(pol, rem, x0, dom)
    end
    println("--------__")
    return vTM
end
