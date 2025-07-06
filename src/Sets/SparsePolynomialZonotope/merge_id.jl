"""
    merge_id(id1::AbstractVector{Int}, id2::AbstractVector{Int}, 
             E₁::AbstractMatrix{N}, E₂::AbstractMatrix{N}) where {N}

Align two identifier vectors and their corresponding exponent matrices in compatible form.

### Input

- `id1::AbstractVector{Int}` -- identifiers corresponding to the rows of `E₁`
- `id2::AbstractVector{Int}` -- identifiers corresponding to the rows of `E₂`
- `E₁::AbstractMatrix{N}` -- first exponent matrix of size `(p₁ × h₁)`
- `E₂::AbstractMatrix{N}` -- second exponent matrix of size `(p₂ × h₂)`

### Output

- `Ē₁::Matrix{N}`: Aligned version of `E₁` 
- `Ē₂::Matrix{N}`: Aligned version of `E₂` 
- `idx::Vector{Int}`: Merged identifier vector, containing all elements of `id1` and any new ones from `id2`.

### Algorithm

This method implements [KochdumperA21; Proposition 1](@citet).

### Example

```jldoctest
julia> import LazySets.SparsePolynomialZonotopeModule: merge_id

julia> id1 = [1, 2];

julia> E₁ = [1 2; 1 0];

julia> id2 = [2, 3];

julia> E₂ = [1 0 1; 3 2 0];

julia> Ē₁, Ē₂, idx = merge_id(id1, id2, E₁, E₂)
([1 2; 1 0; 0 0], [0 0 0; 1 0 1; 3 2 0], [1, 2, 3])
```
"""
function merge_id(id1::AbstractVector{Int}, id2::AbstractVector{Int}, E₁::AbstractMatrix{N},
                  E₂::AbstractMatrix{N}) where {N<:Integer}
    p₁, h₁ = size(E₁)
    p₂, h₂ = size(E₂)
    if !(length(id1) == p₁ && length(id2) == p₂)
        throw(ArgumentError("incompatible dimensions"))
    end

    if id1 == id2
        return E₁, E₂, id1
    end

    # Indices of elements of id2 which do not belong to id1.
    K = findall(∈(setdiff(id2, id1)), id2)
    k = length(K)

    idx = vcat(id1, id2[K])
    Ē₁ = vcat(E₁, zeros(N, k, h₁))

    Ē₂ = zeros(N, p₁ + k, h₂)
    @inbounds for i in 1:(p₁ + k)
        j = findfirst(==(idx[i]), id2)
        if !isnothing(j)
            Ē₂[i, :] = E₂[j, :]
        end
    end
    return Ē₁, Ē₂, idx
end
