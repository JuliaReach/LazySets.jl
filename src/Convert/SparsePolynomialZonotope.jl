"""
    convert(::Type{SparsePolynomialZonotope}, Z::AbstractZonotope{N}) where {N}

Convert a zonotope to sparse polynomial zonotope.

### Input

- `SparsePolynomialZonotope` -- target type
- `Z`                        -- zonotopic set

### Output

A sparse polynomial zonotope.

### Algorithm

The method implements [Kochdumper21a; Proposition 3.1.19](@citet).
"""
function convert(::Type{SparsePolynomialZonotope}, Z::AbstractZonotope{N}) where {N}
    c = center(Z)
    G = genmat(Z)
    p = ngens(Z)
    E = Matrix(1 * I, p, p)
    GI = zeros(N, dim(Z), 0)
    return SparsePolynomialZonotope(c, G, GI, E)
end

"""
    convert(::Type{SparsePolynomialZonotope}, SSPZ::SimpleSparsePolynomialZonotope)

Convert a simple sparse polynomial zonotope to a sparse polynomial zonotope.

### Input

- `SparsePolynomialZonotope` -- target type
- `SSPZ`                     -- simple sparse polynomial zonotope

### Output

A sparse polynomial zonotope.
"""
function convert(::Type{SparsePolynomialZonotope},
                 SSPZ::SimpleSparsePolynomialZonotope{N}) where {N}
    c = center(SSPZ)
    G = genmat(SSPZ)
    E = expmat(SSPZ)
    GI = Matrix{N}(undef, dim(SSPZ), 0)
    return SparsePolynomialZonotope(c, G, GI, E)
end

function load_TaylorModels_convert_SparsePolynomialZonotope()
    return quote
        using .TaylorModels: TaylorModelN, polynomial, remainder, constant_term
        using .TaylorModels.TaylorSeries: coeff_table

        # implements Proposition 3.1.12 in thesis
        function convert(::Type{SparsePolynomialZonotope},
                         vTM::Vector{<:TaylorModelN{r,N}}) where {r,N}
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
                c[i] = IA.mid(rem) + constant_term(pol)
                GI_diag[i] = IA.radius(rem)
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
    end
end  # quote / load_TaylorModels_convert_SparsePolynomialZonotope
