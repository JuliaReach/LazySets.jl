"""
    MatrixZonotope{N, MN<:AbstractMatrix{N}}(A0::MN, Ai::Vector{MN},
                    idx::Vector{Int}=collect(1:length(Aᵢ))) <: AbstractMatrixZonotope{N}

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
Matrix zonotopes were introduced in [AlthoffGK11](@citet).

### Examples

```jldoctest
julia> A0 = [2.0 1.0; -1.0 0.0];

julia> Ai = [[1.0 -1.0; 0.0 -1.0], [0.0 2.0; -1.0 1.0]];

julia> idx = [1, 3];

julia> MZ = MatrixZonotope(A0, Ai, idx)
MatrixZonotope{Float64, Matrix{Float64}}([2.0 1.0; -1.0 0.0], [[1.0 -1.0; 0.0 -1.0], [0.0 2.0; -1.0 1.0]], [1, 3])
```
"""
struct MatrixZonotope{N,MN<:AbstractMatrix{N}} <: AbstractMatrixZonotope{N}
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

function eltype(::Type{<:MatrixZonotope{N}}) where {N}
    return N
end

function size(MZ::MatrixZonotope)
    return size(center(MZ))
end

function size(MZ::MatrixZonotope, d::Int)
    return size(center(MZ), d)
end

"""
    transpose(MZ::MatrixZonotope)

Return the transpose of a matrix zonotope.

### Notes

The transpose of a matrix zonotope is defined as:

```math
    \\mathcal{A}ᵀ = \\braket{(A^{(0)})ᵀ,(A^{(1)})ᵀ, \\dots, (A^{(p)})ᵀ }
```
"""
function transpose(MZ::MatrixZonotope{N}) where {N}
    Ct = transpose(center(MZ))
    Gts = map(transpose, generators(MZ))
    return MatrixZonotope(Ct, Gts)
end

function copy(MZ::MatrixZonotope)
    return MatrixZonotope(copy(MZ.A0),
                          [copy(Aij) for Aij in MZ.Ai],
                          copy(MZ.idx))
end
