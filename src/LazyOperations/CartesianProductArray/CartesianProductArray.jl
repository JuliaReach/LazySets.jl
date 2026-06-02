import Base: *

"""
    CartesianProductArray{N, S<:LazySet{N}} <: LazySet{N}

Type that represents the Cartesian product of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

The Cartesian product preserves convexity: if the set arguments are convex, then
their Cartesian product is convex as well.
"""
struct CartesianProductArray{N,S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

# constructor for an empty product with optional size hint and numeric type
function CartesianProductArray(n::Int=0, N::Type=Float64)
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return CartesianProductArray(arr)
end

"""
    CartesianProduct!(X, Y)

Convenience function to compute the lazy Cartesian product and modify
`CartesianProductArray`s in-place.
"""
function CartesianProduct! end

"""
```
    *(X::LazySet, Xs::LazySet...)
    *(Xs::Vector{<:LazySet})
```

Alias for the n-ary Cartesian product.
"""
*(X::LazySet, Xs::LazySet...) = CartesianProductArray(vcat(X, Xs...))
*(X::LazySet) = X
*(Xs::Vector{<:LazySet}) = CartesianProductArray(Xs)

"""
    ×(X::LazySet, Xs::LazySet...)
    ×(Xs::Vector{<:LazySet})

Alias for the n-ary Cartesian product.

### Notes

The function symbol can be typed via `\\times<tab>`.
"""
×(X::LazySet, Xs::LazySet...) = *(X, Xs...)
×(Xs::Vector{<:LazySet}) = *(Xs)

# add functions connecting CartesianProduct and CartesianProductArray
@declare_array_version(CartesianProduct, CartesianProductArray)

"""
    array(cpa::CartesianProductArray)

Return the array of a Cartesian product of a finite number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

The array of a Cartesian product of a finite number of sets.
"""
function array(cpa::CartesianProductArray)
    return cpa.array
end

"""
    same_block_structure(x::AbstractVector{S1}, y::AbstractVector{S2}
                        ) where {S1<:LazySet, S2<:LazySet}

Check whether two vectors of sets have the same block structure, i.e., the
``i``-th entry in the vectors have the same dimension.

### Input

- `x` -- first vector
- `y` -- second vector

### Output

`true` iff the vectors have the same block structure.
"""
function same_block_structure(x::AbstractVector{S1},
                              y::AbstractVector{S2}) where {S1<:LazySet,S2<:LazySet}
    if length(x) != length(y)
        return false
    end
    for i in eachindex(x)
        if dim(x[i]) != dim(y[i])
            return false
        end
    end
    return true
end

"""
    block_structure(cpa::CartesianProductArray)

Compute an array containing the dimension ranges of each block of a Cartesian
product of a finite number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

A vector of ranges.

### Example

```jldoctest
julia> using LazySets: block_structure

julia> cpa = CartesianProductArray([BallInf(zeros(n), 1.0) for n in [3, 1, 2]]);

julia> block_structure(cpa)
3-element Vector{UnitRange{Int64}}:
 1:3
 4:4
 5:6
```
"""
function block_structure(cpa::CartesianProductArray)
    result = Vector{UnitRange{Int}}(undef, length(array(cpa)))
    start_index = 1
    @inbounds for (i, bi) in enumerate(array(cpa))
        end_index = start_index + dim(bi) - 1
        result[i] = start_index:end_index
        start_index = end_index + 1
    end
    return result
end

"""
    block_to_dimension_indices(cpa::CartesianProductArray{N},
                               vars::Vector{Int}) where {N}

Compute a vector mapping block index `i` to tuple `(f, l)` such that either
`f = l = -1` or `f` is the first dimension index and `l` is the last dimension
index of the `i`-th block, depending on whether one of the block's dimension
indices is specified in `vars`.

### Input

- `cpa`  -- Cartesian product of a finite number of sets
- `vars` -- list containing the variables of interest, sorted in ascending order

### Output

(i) A vector of pairs, where each pair corresponds to the range of dimensions
in the i-th block.

(ii) The number of constrained blocks.

### Example

```jldoctest
julia> using LazySets: block_to_dimension_indices

julia> cpa = CartesianProductArray([BallInf(zeros(n), 1.0) for n in [1, 3, 2, 3]]);

julia> m, k = block_to_dimension_indices(cpa, [2, 4, 8]);

julia> m
4-element Vector{Tuple{Int64, Int64}}:
 (-1, -1)
 (2, 4)
 (-1, -1)
 (7, 9)

julia> k
2
```

The vector `m` represents the mapping "second block from dimension 2 to
dimension 4, fourth block from dimension 7 to dimension 9."
These blocks contain the dimensions specified in `vars=[2, 4, 8]`.
The number of constrained blocks is `k` = 2 (2nd and 4th blocks).
"""
function block_to_dimension_indices(cpa::CartesianProductArray{N},
                                    vars::Vector{Int}) where {N}
    ranges = fill((-1, -1), length(array(cpa)))
    constrained_blocks = 0
    start_index, end_index = 1, 0
    v_i = 1
    @inbounds for i in eachindex(cpa.array)
        end_index += dim(cpa.array[i])
        if v_i <= length(vars) && vars[v_i] <= end_index
            ranges[i] = (start_index, end_index)
            constrained_blocks += 1
            while v_i <= length(vars) && vars[v_i] <= end_index
                v_i += 1
            end
        end
        if v_i > length(vars)
            break
        end
        start_index = end_index + 1
    end
    return ranges, constrained_blocks
end

# method for all variables
function block_to_dimension_indices(cpa::CartesianProductArray{N}) where {N}
    ranges = Vector{Tuple{Int,Int}}(undef, length(cpa.array))

    start_index, end_index = 1, 0
    @inbounds for i in eachindex(cpa.array)
        end_index += dim(cpa.array[i])
        ranges[i] = (start_index, end_index)
        start_index = end_index + 1
    end
    constrained_blocks = length(cpa.array)
    return ranges, constrained_blocks
end

"""
    substitute_blocks(low_dim_cpa::CartesianProductArray,
                      orig_cpa::CartesianProductArray,
                      blocks::Vector{Tuple{Int, Int}})

Return a Cartesian product of a finite number of sets (CPA) obtained by merging
an original CPA with a low-dimensional CPA, which represents the updated subset
of variables in the specified blocks.

### Input

- `low_dim_cpa` -- low-dimensional Cartesian product of a finite number of sets
- `orig_cpa`    -- original high-dimensional Cartesian product of a finite
                   number of sets
- `blocks`      -- index of the first variable in each block of `orig_cpa`

### Output

The merged Cartesian product.
"""
function substitute_blocks(low_dim_cpa::CartesianProductArray,
                           orig_cpa::CartesianProductArray,
                           blocks::Vector{Tuple{Int,Int}})
    N = promote_type(eltype(low_dim_cpa), eltype(orig_cpa))
    array = Vector{LazySet{N}}(undef, length(orig_cpa.array))
    index = 1
    for bi in eachindex(orig_cpa.array)
        start_ind, _ = blocks[bi]
        if start_ind == -1
            array[bi] = orig_cpa.array[bi]
        else
            array[bi] = low_dim_cpa.array[index]
            index += 1
        end
    end
    return CartesianProductArray(array)
end

include("center.jl")
include("concretize.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("ispolyhedral.jl")
include("ispolyhedraltype.jl")
include("vertices_list.jl")
include("volume.jl")
include("in.jl")
include("linear_map.jl")
include("project.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
