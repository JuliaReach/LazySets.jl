module LazySetsIntervalMatricesExt

using IntervalArithmetic: Interval, interval
using IntervalMatrices: AbstractIntervalMatrix, IntervalMatrix, midpoint_radius,
                        _exp_remainder
using LazySets: AbstractZonotope, LinearMap, center, dim, generators,
                minkowski_sum, ngens, scale
using LazySets.ZonotopeModule: Zonotope
using LazySets.MatrixZonotopeModule: AbstractMatrixZonotope, MatrixZonotope,
                                     MatrixZonotopeExp, MatrixZonotopeProduct
using LinearAlgebra: I, dot
using ReachabilityBase.Commutative: @commutative
import Base: convert
import LazySets: minkowski_sum, remove_redundant_generators
import LazySets.Approximations: overapproximate, taylor_expmap_remainder

"""
    taylor_expmap_remainder(Z::AbstractZonotope{N}, matnorm::Real, Int) where {N}

Overapproximate the Lagrange remainder term of the k-th order truncated Taylor expansion
of the exponential map of a matrix zonotope applied to a zonotopic set.

### Input

- `P` -- a zonotopic set
- `matnorm` -- an upper bound on the norm of the matrix zonotope
- `k` -- the order of the Taylor expansion

### Output

A zonotope over-approximating the remainder term of the Taylor expansion.
"""
function taylor_expmap_remainder(Z::AbstractZonotope{N}, matnorm::Real, k::Int) where {N}
    n = dim(Z)
    ϵ = matnorm / (k + 2)

    E = IntervalMatrix(fill(interval(N(-1.0), N(1.0)), n, n))
    E *= matnorm^(k + 1) / (factorial(k + 1) * (1 - ϵ))

    res = overapproximate(E * Z, Zonotope)
    return res
end

# TODO temporary patch for IntervalArithmetic#317 (fixed and redundant in v1)
function convert(::Type{Interval{T}}, x::Interval{T}) where {T<:Real}
    return x
end

# TODO avoid `convert`
@commutative function minkowski_sum(MZ::MatrixZonotope, IM::IntervalMatrix)
    return minkowski_sum(MZ, convert(MatrixZonotope, IM))
end

"""
    overapproximate(lm::LinearMap{N, <:AbstractZonotope, NM,
                                    <:AbstractIntervalMatrix{NM}},
                    ::Type{<:Zonotope}) where {N, NM}

Overapproximate an interval-matrix linear map of a zonotopic set by a zonotope.

### Input

- `lm`       -- interval-matrix linear map of a zonotopic set
- `Zonotope` -- target set type

### Output

A zonotope overapproximating the linear map.

### Algorithm

This implementation uses the method proposed in [AlthoffSB07](@citet).

Given an interval matrix ``M = \\tilde{M} + ⟨-\\hat{M},\\hat{M}⟩`` (split into a
conventional matrix and a symmetric interval matrix) and a zonotope
``⟨c, g_1, …, g_m⟩``, we compute the resulting zonotope
``⟨\\tilde{M}c, \\tilde{M}g_1, …, \\tilde{M}g_m, v_1, …, v_n⟩`` where the
``v_j``, ``j = 1, …, n``, are defined as

```math
    v_j = \\begin{cases} 0 & i ≠ j \\\\
            \\hat{M}_j (|c| + ∑_{k=1}^m |g_k|) & i = j. \\end{cases}
```
"""
function overapproximate(lm::LinearMap{N,<:AbstractZonotope,NM,<:AbstractIntervalMatrix{NM}},
                         ::Type{<:Zonotope}) where {N,NM}
    Mc, Ms = midpoint_radius(lm.M)
    Z = lm.X
    c = Mc * center(Z)
    n = dim(lm)
    nG = ngens(Z)
    G = zeros(N, n, nG + n)
    vector_sum = abs.(center(Z))
    @inbounds for (j, g) in enumerate(generators(Z))
        G[:, j] = Mc * g
        vector_sum += abs.(g)
    end
    @inbounds for i in 1:n
        row = @view Ms[i, :]
        G[i, i + nG] = dot(row, vector_sum)
    end
    return Zonotope(c, G)
end

####################
# matrix zonotopes #
####################

"""
    convert(::Type{MatrixZonotope}, IM::IntervalMatrix)

Convert an interval matrix to a matrix zonotope

### Input

- `MatrixZonotope` -- target type
- `IM` -- an interval matrix

### Output

A matrix zonotope with one generator

### Example

```jldoctest
julia> using LazySets, IntervalMatrices

julia> IM = IntervalMatrix([interval(-1.1, -0.9) interval(-4.1, -3.9);
            interval(3.9, 4.1) interval(-1.1, -0.9)])
2×2 IntervalMatrix{Float64, IntervalArithmetic.Interval{Float64}, Matrix{IntervalArithmetic.Interval{Float64}}}:
 [-1.10001, -0.9]  [-4.1, -3.9]
  [3.9, 4.1]       [-1.10001, -0.9]

julia> MZ = convert(MatrixZonotope, IM)
MatrixZonotope{Float64, Matrix{Float64}}([-1.0 -4.0; 4.0 -1.0], [[0.10000000000000009 0.0; 0.0 0.0], [0.0 0.0; 0.10000000000000009 0.0], [0.0 0.10000000000000009; 0.0 0.0], [0.0 0.0; 0.0 0.10000000000000009]], [1, 2, 3, 4])
```
"""
function convert(::Type{MatrixZonotope}, IM::IntervalMatrix{N}) where {N}
    m, n = size(IM)
    center, halfIM = midpoint_radius(IM)

    gens = Vector{Matrix{N}}()
    sizehint!(gens, m * n)

    @inbounds for (I, val) in enumerate(halfIM)
        iszero(val) && continue
        G = zeros(N, m, n)
        G[I] = val
        push!(gens, G)
    end

    return MatrixZonotope(center, gens)
end

"""
    overapproximate(expA::MatrixZonotopeExp{N,T}, ::Type{<:MatrixZonotope},
                            k::Int=2) where {N,T<:AbstractMatrixZonotope{N}}

Overapproximate the matrix zonotope exponential ``exp(\\mathcal{A})``

### Input

- `expA` -- `MatrixZonotopeExp`
- `MatrixZonotope` -- target type
- `k` -- (default: `2`) the order of the Taylor expansion
- `tol` -- (default: `1e-9`) tolerance used when pruning generators after the overapproximation

### Output

A matrix zonotope overapproximating the matrix zonotope exponential

### Algorithm

The expansion

```math
exp(\\mathcal{A}) ⊆ \\sum_i^k \\frac{\\mathcal{A}^i}{i!} + E_k
```

is computed by overapproximating the matrix zonotope powers ``A^i``
for ``i=0, …, k``.
The remainder term ``E_k`` is computed through interval arithmetic
following [AlthoffKS11; Proposition 4.1](@citet).
"""
function overapproximate(expA::MatrixZonotopeExp{N,T},
                            ::Type{<:MatrixZonotope},
                            k::Int=2; tol::Real=1e-9) where {N,T<:AbstractMatrixZonotope{N}}

    # overapproximate the product MZP = A*B*... ---
    MZP = MatrixZonotopeProduct(expA.M)
    MZ = overapproximate(MZP, MatrixZonotope)

    Id = MatrixZonotope(Matrix{N}(I, size(MZ)), Matrix{N}[])
    powers = Vector{typeof(MZ)}(undef, k)
    powers[1] = MZ
    @inbounds for i in 2:k
        term = overapproximate(MZ * powers[i - 1], MatrixZonotope)
        powers[i] = scale(1 / i, term)
    end
    W = reduce(minkowski_sum, powers)
    W = minkowski_sum(W, Id)

    IM = overapproximate(MZ, IntervalMatrix)
    E = _exp_remainder(IM, N(1), k)
    res = minkowski_sum(W, E)

    return remove_redundant_generators(res; tol=tol)
end

"""
    overapproximate(A::MatrixZonotope, ::Type{<:IntervalMatrix})

Overapproximate a matrix zonotope with an interval matrix.

### Input

- `A` -- a matrix zonotope
- `IntervalMatrix` -- target type

### Output

An interval matrix overapproximating the matrix zonotope.
"""
function overapproximate(A::MatrixZonotope, ::Type{<:IntervalMatrix})
    A_abs = sum([abs.(Ai) for Ai in generators(A)])
    A₊ = center(A) + A_abs
    A₋ = center(A) - A_abs
    return IntervalMatrix(A₋, A₊)
end

end  # module
