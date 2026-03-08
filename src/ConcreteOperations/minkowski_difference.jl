"""
    minkowski_difference(P::LazySet, Q::LazySet)

Concrete Minkowski difference (geometric difference) of a polytopic set and a
compact set.

### Input

- `P` -- polytopic set
- `Q` -- compact set that is subtracted from `P`

### Output

An `HPolytope` that corresponds to the Minkowski difference of `P` minus `Q` if
`P` is bounded, and an `HPolyhedron` if `P` is unbounded.

### Notes

This implementation requires that the set `P` is polyhedral and that the set `Q`
is bounded.

### Algorithm

This method implements [KolmanovskyG98; Theorem 2.3](@citet):

Suppose ``P`` is a polyhedron
```math
P = \\{z ∈ ℝ^n: sᵢᵀz ≤ rᵢ,~i = 1, …, k\\}.
```
where ``sᵢ ∈ ℝ^n, sᵢ ≠ 0``, and ``rᵢ ∈ ℝ``.
Assume ``ρ(sᵢ,Q)`` is defined for ``i = 1, …, k``.
Then the Minkowski difference is

```math
\\{z ∈ ℝ^n: sᵢᵀz ≤ rᵢ - ρ(sᵢ,Q),~i = 1, …, k\\}.
```

While the algorithm applies the support function to `Q`, we have that
``P ⊖ Q = P ⊖ \\text{CH}(Q)`` whenever `P` is convex, where CH denotes the
convex hull. Hence, if `Q` is not convex by type information, we wrap it in a
lazy `ConvexHull`.
"""
@validate function minkowski_difference(P::LazySet, Q::LazySet)
    @assert ispolyhedral(P) "this implementation requires that the first argument " *
                            "is polyhedral; try overapproximating with an `HPolyhedron`"
    @assert isbounded(Q) "this implementation requires that the second " *
                         "argument is bounded, but it is not"

    if !isconvextype(typeof(Q))
        Q = ConvexHull(Q)
    end

    A, b = tosimplehrep(P)
    N = promote_type(eltype(P), eltype(Q))
    g_PminusQ = Vector{N}(undef, length(b))
    @inbounds for i in eachindex(b)
        g_PminusQ[i] = b[i] - ρ(A[i, :], Q)
    end
    if isbounded(P)
        return HPolytope(A, g_PminusQ)
    else
        return HPolyhedron(A, g_PminusQ)
    end
end

for T in (:LazySet, :AbstractZonotope, :AbstractHyperrectangle)
    # Minkowski difference with singleton is a translation
    @eval @validate function minkowski_difference(X::($T), S::AbstractSingleton)
        return translate(X, -element(S))
    end

    # Minkowski difference with ZeroSet is the identity
    @eval @validate function minkowski_difference(X::($T), Z::ZeroSet)
        return X
    end
end

"""
# Extended help

    minkowski_difference(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle)

### Output

A `Hyperrectangle`, or an `EmptySet` if the difference is empty.
"""
@validate function minkowski_difference(H1::AbstractHyperrectangle,
                                        H2::AbstractHyperrectangle)
    n = dim(H1)
    N = promote_type(eltype(H1), eltype(H2))
    r = Vector{N}(undef, n)
    for i in 1:n
        r[i] = radius_hyperrectangle(H1, i) - radius_hyperrectangle(H2, i)
        if r[i] < zero(N)
            return EmptySet{N}(n)
        end
    end
    return Hyperrectangle(center(H1) - center(H2), r)
end

"""
# Extended help

    minkowski_difference(Z1::AbstractZonotope, Z2::AbstractZonotope)

### Output

An `HPolytope` for dimensions greater than two, and a `Zonotope` or `EmptySet`
for lower dimensions.

### Algorithm

For one-dimensional sets, this method implements a simple algorithm for intervals.
For two-dimensional sets, this method implements [Althoff15; Proposition 6](@citet).
For higher-dimensional sets, this method implements [Althoff15; Theorem 3](@citet).
"""
@validate function minkowski_difference(Z1::AbstractZonotope, Z2::AbstractZonotope)
    n = dim(Z1)
    if n == 1
        return _minkowski_difference_1d(Z1, Z2)
    elseif n == 2
        return _minkowski_difference_2d(Z1, Z2)
    else
        return _minkowski_difference_nd(Z1, Z2)
    end
end

function _minkowski_difference_1d(Z1::AbstractZonotope, Z2::AbstractZonotope)
    N = promote_type(eltype(I1), eltype(I2))
    l = low(I1, 1) - low(I2, 1)
    if h < l
        return EmptySet{N}(1)
    end
    c = center(Z1, 1) - center(Z2, 1)
    return Zonotope([c], hcat([(c - l) / 2]))
end

# See [Althoff15; Proposition 6](@citet).
function _minkowski_difference_2d(Zm::AbstractZonotope, Zs::AbstractZonotope)
    N = promote_type(eltype(Zm), eltype(Zs))
    Gm = genmat(Zm)
    n, p = size(Gm)
    sii = StrictlyIncreasingIndices(p, n - 1)
    C⁺ = Matrix{N}(undef, length(sii), n)
    for (i, columns) in enumerate(sii)
        c⁺ = cross_product(view(Gm, :, columns))
        normalize!(c⁺, 2)
        C⁺[i, :] .= c⁺
    end
    Δd = sum(g -> abs.(C⁺ * g), generators(Zm))
    Δdtrans = sum(g -> abs.(C⁺ * g), generators(Zs))
    # stretching factors
    μ = inv(abs.(C⁺ * Gm)) * (Δd .- Δdtrans)
    # if a factor is zero, the generator will be zero and can be removed
    # if a factor is negative, the resulting set is empty
    zero_gens = 0
    for v in μ
        if isapproxzero(v)
            zero_gens += 1
        elseif v < zero(N)
            return EmptySet{N}(2)
        end
    end
    # stretch generators g_j of Gm by μ[j] (for all j) and skip zero cases
    G = Matrix{N}(undef, n, p - zero_gens)
    col = 1
    for j in 1:p
        if !isapproxzero(μ[j])
            G[:, col] .= Gm[:, j] * μ[j]
            col += 1
        end
    end
    c = center(Zm) .- center(Zs)
    return Zonotope(c, G)
end

# See [Althoff15; Theorem 3](@citet).
function _minkowski_difference_nd(Z1::AbstractZonotope, Z2::AbstractZonotope)
    Gm = genmat(Z1)
    n, p = size(Gm)

    N = promote_type(eltype(Z1), eltype(Z2))
    cm, Gmᵀ = center(Z1), transpose(Gm)
    cs, Gsᵀ = center(Z2), transpose(genmat(Z2))
    Δc = cm - cs

    constraints = Vector{HalfSpace{N,Vector{N}}}()
    for columns in StrictlyIncreasingIndices(p, n - 1)
        c⁺ = cross_product(view(Gm, :, columns))
        iszero(c⁺) && continue
        normalize!(c⁺, 2)

        Δd = sum(abs, Gmᵀ * c⁺)
        Δdtrans = sum(abs, Gsᵀ * c⁺)

        c⁺Δc = dot(c⁺, Δc)
        ΔΔd = Δd - Δdtrans
        d⁺ = c⁺Δc + ΔΔd
        c⁻ = -c⁺
        d⁻ = -c⁺Δc + ΔΔd

        push!(constraints, HalfSpace(c⁺, d⁺))
        push!(constraints, HalfSpace(c⁻, d⁻))
    end

    return HPolytope(constraints)
end

@validate function minkowski_difference(∅::EmptySet, X::LazySet)
    return _minkowski_difference_emptyset(∅, X)
end

@validate function minkowski_difference(X::LazySet, ∅::EmptySet)
    return _minkowski_difference_emptyset2(X, ∅)
end

@validate function minkowski_difference(U::Universe, X::LazySet)
    return _minkowski_difference_universe(U, X)
end

@validate function minkowski_difference(X::LazySet, U::Universe)
    return _minkowski_difference_universe2(X, U)
end

# ============== #
# disambiguation #
# ============== #

for T in (:AbstractSingleton, :ZeroSet, :Universe)
    @eval @validate function minkowski_difference(∅::EmptySet, X::($T))
        return _minkowski_difference_emptyset(∅, X)
    end
end

for T in (:AbstractSingleton, :ZeroSet)
    @eval @validate function minkowski_difference(U::Universe, X::$T)
        return _minkowski_difference_universe(U, X)
    end
end

@validate function minkowski_difference(U::Universe, ∅::EmptySet)
    return _minkowski_difference_emptyset2(U, ∅)
end
