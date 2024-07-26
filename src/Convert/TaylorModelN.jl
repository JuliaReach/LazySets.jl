function load_TaylorModels_convert_TaylorModelN()
    return quote
        using .TaylorModels: TaylorModelN
        using .TaylorModels.TaylorSeries: coeff_table, set_variables, HomogeneousPolynomial,
                                          in_base, pos_table, TaylorN

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
    end
end  # quote / load_TaylorModels_convert_TaylorModelN
