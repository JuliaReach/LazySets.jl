import Base: *

export CartesianProductArray,
       array,
       same_block_structure

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

isoperationtype(::Type{<:CartesianProductArray}) = true
isconvextype(::Type{CartesianProductArray{N,S}}) where {N,S} = isconvextype(S)

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
    dim(cpa::CartesianProductArray)

Return the dimension of a Cartesian product of a finite number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

The ambient dimension of the Cartesian product of a finite number of sets, or
`0` if there is no set in the array.
"""
function dim(cpa::CartesianProductArray)
    return length(cpa.array) == 0 ? 0 : sum(dim(Xi) for Xi in cpa.array)
end

"""
    σ(d::AbstractVector, cpa::CartesianProductArray)

Compute a support vector of a Cartesian product of a finite number of sets.

### Input

- `d`   -- direction
- `cpa` -- Cartesian product of a finite number of sets

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the product sets.
"""
function σ(d::AbstractVector, cpa::CartesianProductArray)
    svec = similar(d)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        svec[i0:i1] = σ(d[i0:i1], Xi)
        i0 = i1 + 1
    end
    return svec
end

# faster version for sparse vectors
function σ(d::AbstractSparseVector, cpa::CartesianProductArray)
    # idea: We walk through the blocks of `cpa` (i.e., the sets `Xi`) and search
    # for corresponding non-zero entries in `d` (stored in `indices`).
    # `next_idx` is the next index of `indices` such that
    # `next_dim = indices[next_idx]` lies in the next block to consider
    # (potentially skipping some blocks).
    svec = similar(d)
    indices, _ = SparseArrays.findnz(d)
    if isempty(indices)
        # direction is the zero vector
        return an_element(cpa)
    end
    next_idx = 1
    next_dim = indices[next_idx]
    m = length(indices)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        if next_dim <= i1
            # there is a non-zero entry in this block
            svec[i0:i1] = σ(d[i0:i1], Xi)

            # find next index outside the current block
            next_idx += 1
            while next_idx <= m && indices[next_idx] <= i1
                next_idx += 1
            end
            if next_idx <= m
                next_dim = indices[next_idx]
            end
        else
            svec[i0:i1] = an_element(Xi)
        end
        i0 = i1 + 1
    end
    return svec
end

# faster version for single-entry vectors
function σ(d::SingleEntryVector, cpa::CartesianProductArray)
    svec = similar(d)
    i0 = 1
    idx = d.i
    for Xi in cpa.array
        ni = dim(Xi)
        i1 = i0 + ni - 1
        if i0 <= idx && idx <= i1
            svec[i0:i1] = σ(SingleEntryVector(d.i - i0 + 1, ni, d.v), Xi)
        else
            svec[i0:i1] = an_element(Xi)
        end
        i0 = i1 + 1
    end
    return svec
end

"""
    ρ(d::AbstractVector, cpa::CartesianProductArray)

Evaluate the support function of a Cartesian product of a finite number of sets.

### Input

- `d`   -- direction
- `cpa` -- Cartesian product of a finite number of sets

### Output

The evaluation of the support function in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
function ρ(d::AbstractVector, cpa::CartesianProductArray)
    N = promote_type(eltype(d), eltype(cpa))
    sfun = zero(N)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        sfun += ρ(d[i0:i1], Xi)
        i0 = i1 + 1
    end
    return sfun
end

# faster version for sparse vectors
function ρ(d::AbstractSparseVector, cpa::CartesianProductArray)
    N = promote_type(eltype(d), eltype(cpa))
    # idea: see the σ method for AbstractSparseVector
    sfun = zero(N)
    indices, _ = SparseArrays.findnz(d)
    if isempty(indices)
        # direction is the zero vector
        return sfun
    end
    next_idx = 1
    next_dim = indices[next_idx]
    m = length(indices)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        if next_dim <= i1
            # there is a non-zero entry in this block
            sfun += ρ(d[i0:i1], Xi)

            # find next index outside the current block
            next_idx += 1
            while next_idx <= m && indices[next_idx] <= i1
                next_idx += 1
            end
            if next_idx > m
                # no more non-zero entries
                break
            end
            next_dim = indices[next_idx]
        end
        i0 = i1 + 1
    end
    return sfun
end

# faster version for single-entry vectors
function ρ(d::SingleEntryVector, cpa::CartesianProductArray)
    i0 = 1
    idx = d.i
    for Xi in cpa.array
        ni = dim(Xi)
        i1 = i0 + ni - 1
        if i0 <= idx && idx <= i1
            return ρ(SingleEntryVector(d.i - i0 + 1, ni, d.v), Xi)
        end
        i0 = i1 + 1
    end
    return sfun
end

"""
    isbounded(cpa::CartesianProductArray)

Check whether a Cartesian product of a finite number of sets is bounded.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cpa::CartesianProductArray)
    return all(isbounded, cpa.array)
end

function isboundedtype(::Type{<:CartesianProductArray{N,S}}) where {N,S}
    return isboundedtype(S)
end

function ispolyhedral(cpa::CartesianProductArray)
    return all(ispolyhedral, array(cpa))
end

"""
    ∈(x::AbstractVector, cpa::CartesianProductArray)

Check whether a given point is contained in a Cartesian product of a finite
number of sets.

### Input

- `x`   -- point/vector
- `cpa` -- Cartesian product of a finite number of sets

### Output

`true` iff ``x ∈ \\text{cpa}``.
"""
function ∈(x::AbstractVector, cpa::CartesianProductArray)
    @assert length(x) == dim(cpa)

    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        if x[i0:i1] ∉ Xi
            return false
        end
        i0 = i1 + 1
    end
    return true
end

"""
    isempty(cpa::CartesianProductArray)

Check whether a Cartesian product of a finite number of sets is empty.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

`true` iff any of the sub-blocks is empty.
"""
function isempty(cpa::CartesianProductArray)
    return any(isempty, array(cpa))
end

"""
    center(cpa::CartesianProductArray)

Compute the center of a Cartesian product of a finite number of
centrally-symmetric sets.

### Input

- `cpa` -- Cartesian product of a finite number of centrally-symmetric sets

### Output

The center of the Cartesian product of a finite number of sets.
"""
function center(cpa::CartesianProductArray)
    return reduce(vcat, center(X) for X in cpa)
end

"""
    constraints_list(cpa::CartesianProductArray)

Compute a list of constraints of a (polyhedral) Cartesian product of a finite
number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

A list of constraints.
"""
function constraints_list(cpa::CartesianProductArray)
    return _constraints_list_cartesian_product(cpa)
end

function _constraints_list_cartesian_product(cp::Union{CartesianProduct,CartesianProductArray})
    N = eltype(cp)
    clist = Vector{HalfSpace{N,SparseVector{N,Int}}}()
    n = dim(cp)
    sizehint!(clist, n)
    prev_step = 1
    # create high-dimensional constraints list
    for c_low in cp
        c_low_list = constraints_list(c_low)
        if isempty(c_low_list)
            n_low = dim(c_low)
        else
            n_low = dim(c_low_list[1])
            indices = prev_step:(prev_step + n_low - 1)
        end
        for constr in c_low_list
            new_constr = HalfSpace(sparsevec(indices, constr.a, n), constr.b)
            push!(clist, new_constr)
        end
        prev_step += n_low
    end

    return clist
end

"""
    vertices_list(cpa::CartesianProductArray)

Compute a list of vertices of a (polytopic) Cartesian product of a finite
number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

A list of vertices.

### Algorithm

We assume that the underlying sets are polytopic.
Then the high-dimensional set of vertices is just the Cartesian product of the
low-dimensional sets of vertices.
"""
function vertices_list(cpa::CartesianProductArray)
    # collect low-dimensional vertices lists
    vlist_low = [vertices_list(X) for X in cpa]

    # create high-dimensional vertices list
    indices_max = [length(vl) for vl in vlist_low]
    m = prod(indices_max)
    N = eltype(cpa)
    vlist = Vector{Vector{N}}(undef, m)
    indices = ones(Int, length(vlist_low))
    v = zeros(N, dim(cpa))
    dim_start_j = 1
    for vl in vlist_low
        v_low = vl[1]
        v[dim_start_j:(dim_start_j + length(v_low) - 1)] = v_low
        dim_start_j += length(v_low)
    end
    i = 1
    j = 1
    # iterate through all index combinations
    while true
        indices[1] = 0
        while indices[1] < indices_max[1]
            indices[1] += 1
            v_low = vlist_low[1][indices[1]]
            v[1:length(v_low)] = v_low
            vlist[i] = copy(v)
            i += 1
        end
        if i > m
            break
        end
        j = 1
        dim_start_j = 1
        while indices[j] == indices_max[j]
            indices[j] = 1
            v_low = vlist_low[j][1]
            v[dim_start_j:(dim_start_j + length(v_low) - 1)] = v_low
            dim_start_j += length(v_low)
            j += 1
        end
        indices[j] += 1
        v_low = vlist_low[j][indices[j]]
        v[dim_start_j:(dim_start_j + length(v_low) - 1)] = v_low
    end

    return vlist
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
    substitute_blocks(low_dim_cpa::CartesianProductArray{N},
                      orig_cpa::CartesianProductArray{N},
                      blocks::Vector{Tuple{Int, Int}}) where {N}

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
function substitute_blocks(low_dim_cpa::CartesianProductArray{N},
                           orig_cpa::CartesianProductArray{N},
                           blocks::Vector{Tuple{Int,Int}}) where {N}
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

"""
    linear_map(M::AbstractMatrix, cpa::CartesianProductArray)

Concrete linear map of a Cartesian product of a finite number of (polyhedral)
sets.

### Input

- `M`   -- matrix
- `cpa` -- Cartesian product of a finite number of sets

### Output

A polyhedron or polytope.
"""
function linear_map(M::AbstractMatrix, cpa::CartesianProductArray)
    return _linear_map_cartesian_product(M, cpa)
end

function project(cpa::CartesianProductArray, block::AbstractVector{Int};
                 kwargs...)
    target_sets = LazySet[]
    m = length(block)

    # find first set
    i_start = 1
    bi = block[i_start]
    n_sum = 0
    n_sum_old = 0
    @inbounds for (j, Xj) in enumerate(array(cpa))
        nj = dim(Xj)
        n_sum += nj
        if n_sum >= bi
            # found starting point in a set; now find end point
            i_end = m
            for i in (i_start + 1):m
                if block[i] > n_sum
                    i_end = i - 1
                    break
                end
            end

            # project this block
            projected = project(Xj, block[i_start:i_end] .- n_sum_old; kwargs...)
            push!(target_sets, projected)

            if i_end == m
                # last index visited
                break
            end

            # advance indices
            i_start = i_end + 1
            bi = block[i_start]
        end
        n_sum_old = n_sum
    end

    # construct result depending on the number of sets
    if length(target_sets) == 1
        @inbounds return target_sets[1]
    elseif length(target_sets) == 2
        @inbounds return CartesianProduct(target_sets[1], target_sets[2])
    else
        # create a new array for better type information
        return CartesianProductArray([X for X in target_sets])
    end
end

function concretize(cpa::CartesianProductArray)
    a = array(cpa)
    @assert !isempty(a) "an empty Cartesian product is not allowed"
    X = cpa
    @inbounds for (i, Y) in enumerate(a)
        if i == 1
            X = concretize(Y)
        else
            X = cartesian_product(X, concretize(Y))
        end
    end
    return X
end

"""
    volume(cpa::CartesianProductArray)

Compute the volume of a Cartesian product of a finite number of sets.

### Input

- `cpa` -- Cartesian product of a finite number of sets

### Output

The volume.
"""
function volume(cpa::CartesianProductArray)
    return prod(volume, array(cpa))
end

function translate(cpa::CartesianProductArray, x::AbstractVector)
    res = Vector{LazySet}(undef, length(array(cpa)))
    s = 1
    @inbounds for (j, Xj) in enumerate(array(cpa))
        e = s + dim(Xj) - 1
        res[j] = translate(Xj, @view(x[s:e]))
        s = e + 1
    end
    return CartesianProductArray([X for X in res])
end
