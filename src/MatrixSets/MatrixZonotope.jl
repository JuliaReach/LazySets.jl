"""
    MatrixZonotope{N, MN<:AbstractMatrix{N}}(A0::MN, Ai::Vector{MN},
                                              idx::Vector{Int}=collect(1:length(Aᵢ)))

Type that represents a matrix zonotope.

### Fields

- `A0`  -- center of the matrix zonotope
- `Ai`  -- vector of matrices; each matrix is a generator of the matrix zonotope
- `idx` -- identifier vector of positive integers for each factor

### Notes

Mathematically a matrix zonotope is defined as the set of matrices

```math
\\mathcal{A} = \\left\\{A ∈ ℝ^{n×m} : A^{(0)} + ∑_{i=1}^p ξ_i A^{(i)},~~ ξ_i ∈ [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```

It can be written in shorthand notation as ``\\mathcal{A} = \\braket{A^{(0)},A^{(1)}, ..., A^{(p)} }_{MZ}.
Matrix zonotopes were introduced in [HuangLBS2025](@citet).

### Examples

```jldoctest
julia> A0 = [2.0 1.0; -1.0 0.0];

julia> Ai = [[1.0 -1.0; 0.0 -1.0], [0.0 2.0; -1.0 1.0]];

julia> idx = [1, 3];

julia> MZ = MatrixZonotope(A0, Ai, idx)
MatrixZonotope{Float64, Matrix{Float64}}([2.0 1.0; -1.0 0.0], [[1.0 -1.0; 0.0 -1.0], [0.0 2.0; -1.0 1.0]], [1, 3])
```
"""
struct MatrixZonotope{N,MN<:AbstractMatrix{N}}
    A0::MN
    Ai::Vector{MN}
    idx::Vector{Int}

    function MatrixZonotope(A0::MN, Ai::Vector{MN},
                            idx::Vector{Int}=collect(1:length(Ai))) where {N,MN<:AbstractMatrix{N}}
        length(Ai) == length(idx) == length(unique(idx)) ||
            throw(ArgumentError("the number of generator matrices doesn't match the id's length"))
        if length(Ai) > 0
            (all(size(Aij) == size(A0) for Aij in Ai)) ||
                throw(ArgumentError("the size of all generator matrices should match"))
        end
        return new{N,MN}(A0, Ai, idx)
    end
end

Base.eltype(::Type{<:MatrixZonotope{N}}) where {N} = N

"""
    size(MZ::MatrixZonotope)

Return the dimensions of the center matrix of the matrix zonotope `MZ`.  

"""
Base.size(MZ::MatrixZonotope) = size(center(MZ))
Base.size(MZ::MatrixZonotope, d::Int) = size(center(MZ), d)

"""
    transpose(MZ::MatrixZonotope)

Returns the transpose of the matrix zonotope `MZ`.

### Notes

The transpose of a matrix zonotope is defined as:

```math
    \\mathcal{A}^⊺ = \\braket{(A^{(0)})^⊺,(A^{(1)})^⊺, \\dots, (A^{(p)})^⊺ }
```
"""
Base.transpose(MZ::MatrixZonotope{N}) where {N} = begin
    Ct = transpose(center(MZ))
    Gts = map(transpose, generators(MZ))
    MatrixZonotope(Ct, Gts)
end