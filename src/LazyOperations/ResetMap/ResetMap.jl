using ReachabilityBase.Arrays: find_unique_nonzero_entry

"""
    ResetMap{N, S<:LazySet{N}} <: AbstractAffineMap{N, S}

Type that represents a lazy reset map.
A reset map is a special case of an affine map ``A x + b, x ∈ X`` where the
linear map ``A`` is the identity matrix with zero entries in all reset
dimensions, and the translation vector ``b`` is zero in all other dimensions.

### Fields

- `X`      -- set
- `resets` -- resets (a mapping from an index to a new value)

### Notes

The reset map preserves convexity: if `X` is convex, then any reset map of `X`
is convex as well.

### Examples

```jldoctest resetmap
julia> X = BallInf([2.0, 2.0, 2.0], 1.0);

julia> r = Dict(1 => 4.0, 3 => 0.0);

julia> rm = ResetMap(X, r);

```

Here `rm` modifies the set `X` such that `x1` is reset to 4 and `x3` is reset to
0, while `x2` is not modified.
Hence `rm` is equivalent to the set
`VPolytope([[4.0, 1.0, 0.0], [4.0, 3.0, 0.0]])`, i.e., an axis-aligned line
segment embedded in 3D.

The corresponding affine map ``A x + b`` would be:

```math
    \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix} x +
    \begin{pmatrix} 4 & 0 & 0 \end{pmatrix}
```

Use the function `matrix` (resp. `vector`) to create the matrix `A` (resp.
vector `b`) corresponding to a given reset map.

```jldoctest resetmap
julia> matrix(rm)
3×3 LinearAlgebra.Diagonal{Float64, Vector{Float64}}:
 0.0   ⋅    ⋅
  ⋅   1.0   ⋅
  ⋅    ⋅   0.0

julia> vector(rm)
3-element SparseArrays.SparseVector{Float64, Int64} with 1 stored entry:
  [1]  =  4.0
```

The application of a `ResetMap` to a `ZeroSet` or an `EmptySet` is simplified
automatically.

```jldoctest resetmap
julia> ResetMap(ZeroSet(3), r)
Singleton{Float64, SparseArrays.SparseVector{Float64, Int64}}(sparsevec([1], [4.0], 3))

julia> ResetMap(EmptySet(3), r)
∅(3)
```

The (in this case unique) support vector of `rm` in direction `[1, 1, 1]` is:

```jldoctest resetmap
julia> σ(ones(3), rm)
3-element Vector{Float64}:
 4.0
 3.0
 0.0
```
"""
struct ResetMap{N,S<:LazySet{N}} <: AbstractAffineMap{N,S}
    X::S
    resets::Dict{Int,N}
end

# ZeroSet is "almost absorbing" for the reset map because only the translation
# vector remains
function ResetMap(Z::ZeroSet{N}, resets::Dict{Int,N}) where {N}
    return Singleton(_vector_from_dictionary(resets, dim(Z)))
end

# EmptySet is absorbing for ResetMap
function ResetMap(∅::EmptySet{N}, resets::Dict{Int,N}) where {N}
    return ∅
end

"""
    matrix(rm::ResetMap)

Return the ``A`` matrix of the affine map ``A x + b, x ∈ X`` represented by a
reset map.

### Input

- `rm` -- reset map

### Output

The (`Diagonal`) matrix for the affine map ``A x + b, x ∈ X`` represented by the
reset map.

### Algorithm

We construct the identity matrix and set all entries in the reset dimensions to
zero.
"""
function matrix(rm::ResetMap)
    N = eltype(rm)
    n = dim(rm)
    v = ones(N, n)
    for i in keys(rm.resets)
        v[i] = zero(N)
    end
    return Diagonal(v)
end

"""
    vector(rm::ResetMap)

Return the ``b`` vector of the affine map ``A x + b, x ∈ X`` represented by a
reset map.

### Input

- `rm` -- reset map

### Output

The (sparse) vector for the affine map ``A x + b, x ∈ X`` represented by the
reset map.
The vector contains the reset value for all reset dimensions and is zero for all
other dimensions.
"""
function vector(rm::ResetMap)
    return _vector_from_dictionary(rm.resets, dim(rm))
end

function _vector_from_dictionary(dict::Dict{Int,N}, n::Int) where {N}
    b = sparsevec(Int[], N[], n)
    for (i, val) in dict
        b[i] = val
    end
    return b
end

"""
    set(rm::ResetMap)

Return the set wrapped by a reset map.

### Input

- `rm` -- reset map

### Output

The wrapped set.
"""
function set(rm::ResetMap)
    return rm.X
end

include("an_element.jl")
include("concretize.jl")
include("constraints_list.jl")
include("dim.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isoperationtype.jl")
include("support_function.jl")
include("support_vector.jl")
