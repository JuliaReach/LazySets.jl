"""
    struct ZonotopeMD{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, DN<:AbstractVector{N}} <: AbstractZonotope{N}

Type that represents a zonotope of order `k` in normal form.

### Fields

- `center::VN` — the center of the zonotope
- `M::MN` — matrix of general (non-axis-aligned) generators
- `d::DN` — vector of axis-aligned (diagonal) generators

### Notes
A zonotope is of order `k` if it has `n * k` generators in `ℝⁿ`, where `n` is the ambient dimension.

A zonotope of order `k` in *normal form* is defined as the set

```math
Z = \\left\\{ x ∈ ℝ^n : x = c + Mξ + d ⊙ η, ~~ ξ ∈ [-1, 1]^m, ~~ η ∈ [-1, 1]^n \\right\\},
```

where `M ∈ ℝ^{n×m}` is a matrix of general generators with `m = n*(k -1)` and `d ∈ ℝⁿ` is a vector of axis-aligned generators. 
Equivalently, this can be seen as a zonotope with generator matrix `[M  D]`, where `D` is the diagonal matrix 
formed from the vector `d`.
ZonotopeMD can be constructed in two ways: by passing the full generator matrix `[M  D]` in normal form 
or by passing `M` and a vector `d` separately.

### Examples

Constructing a zonotope in normal form from a center, general generator matrix `M`, and diagonal vector `d`:

```jldoctest zonotopemd_label
julia> c = [0.0, 0.0];

julia> M = [1.0 2.0; 3.0 1.0];

julia> d = [0.1, 0.2];

julia> Z = ZonotopeMD(c, M, d)
ZonotopeMD{Float64, Vector{Float64}, Matrix{Float64}, Vector{Float64}}([0.0, 0.0], [1.0 2.0; 3.0 1.0], [0.1, 0.2])

julia> center(Z)
2-element Vector{Float64}:
 0.0
 0.0

julia> genmat(Z)
2×4 SparseArrays.SparseMatrixCSC{Float64, Int64} with 6 stored entries:
 1.0  2.0  0.1   ⋅ 
 3.0  1.0   ⋅   0.1
```

The generator matrix returned by `genmat` is the concatenation `[M D]`, where `D` is the diagonal matrix formed from `d`.
THe resulting matrix is stored as a sparse matrix.
Constructing the same zonotope by passing the full generator matrix `[M D]` directly:

```jldoctest zonotopemd_label
julia> G = [1.0 2.0 0.1 0.0;
            3.0 1.0 0.0 0.2];

julia> Z2 = ZonotopeMD([0.0, 0.0], G)
ZonotopeMD{Float64, Vector{Float64}, Matrix{Float64}, Vector{Float64}}([0.0, 0.0], [1.0 2.0; 3.0 1.0], [0.1, 0.2])

julia> genmat(Z2) == G
true
```
You can also convert back to a standard `Zonotope` if needed:

```jldoctest zonotopemd_label
julia> Zstd = Zonotope(Z)
Zonotope{Float64, Vector{Float64}, SparseArrays.SparseMatrixCSC{Float64, Int64}}([0.0, 0.0], sparse([1, 2, 1, 2, 1, 2], [1, 1, 2, 2, 3, 4], [1.0, 3.0, 2.0, 1.0, 0.1, 0.2], 2, 4))
```

"""
struct ZonotopeMD{N,VN<:AbstractVector{N},MN<:AbstractMatrix{N},DN<:AbstractVector{N}} <:
       AbstractZonotope{N}
    center::VN
    M::MN
    d::DN

    function ZonotopeMD(center::VN, M::MN,
                        d::DN) where {N,
                                      VN<:AbstractVector{N},
                                      MN<:AbstractMatrix{N},
                                      DN<:AbstractVector{N}}
        @assert length(center) == size(M, 1) == length(d) "Dimensions must match"
        return new{N,VN,MN,DN}(center, M, d)
    end
end

# constructor from generator matrix
function ZonotopeMD(center::VN, G::AbstractMatrix{N}) where {N,VN<:AbstractVector{N}}
    n, p = size(G)
    @assert p % n == 0 "The generator matrix must contain a multiple of n columns"
    @assert p >= 2n "Expected at least order 2 zonotope"

    M = G[:, 1:(p - n)]
    D = G[:, (end - n + 1):end]

    @assert isdiag(D) "The last n columns of the generator matrix must be diagonal"
    d = diag(D)
    return ZonotopeMD(center, M, d)
end
