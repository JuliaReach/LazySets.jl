import Base: *, ∈, isempty

export CartesianProduct,
       CartesianProductArray,
       CartesianProduct!,
       array,
       swap,
       same_block_structure

"""
    CartesianProduct{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents a Cartesian product of two convex sets.

### Fields

- `X` -- first convex set
- `Y` -- second convex set

### Notes

The Cartesian product of three elements is obtained recursively.
See also `CartesianProductArray` for an implementation of a Cartesian product of
many sets without recursion, instead using an array.

The `EmptySet` is the absorbing element for `CartesianProduct`.

### Examples

The Cartesian product between two sets `X` and `Y` can be constructed either
using `CartesianProduct(X, Y)` or the short-cut notation `X × Y`:

```jldoctest cartesianproduct_constructor
julia> I1 = Interval(0, 1);

julia> I2 = Interval(2, 4);

julia> I12 = I1 × I2;

julia> typeof(I12)
CartesianProduct{Float64,Interval{Float64,IntervalArithmetic.Interval{Float64}},Interval{Float64,IntervalArithmetic.Interval{Float64}}}
```
A hyperrectangle is the cartesian product of intervals, so we can convert `I12`
exactly to a `Hyperrectangle` type:

```jldoctest cartesianproduct_constructor
julia> convert(Hyperrectangle, I12)
Hyperrectangle{Float64}([0.5, 3.0], [0.5, 1.0])
```
"""
struct CartesianProduct{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2
end

# EmptySet is the absorbing element for CartesianProduct
@absorbing(CartesianProduct, EmptySet)

"""
```
    *(X::LazySet, Y::LazySet)
```

Alias for the binary Cartesian product.
"""
*(X::LazySet, Y::LazySet) = CartesianProduct(X, Y)

"""
    ×

Alias for the binary Cartesian product.
"""
×(X::LazySet, Y::LazySet) = CartesianProduct(X, Y)

"""
    swap(cp::CartesianProduct)

Return a new `CartesianProduct` object with the arguments swapped.

### Input

- `cp` -- Cartesian product of two convex sets

### Output

A new `CartesianProduct` object with the arguments swapped.
"""
function swap(cp::CartesianProduct)
    return CartesianProduct(cp.Y, cp.X)
end

"""
    dim(cp::CartesianProduct)::Int

Return the dimension of a Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

The ambient dimension of the Cartesian product.
"""
function dim(cp::CartesianProduct)::Int
    return dim(cp.X) + dim(cp.Y)
end

"""
    σ(d::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}

Return the support vector of a Cartesian product.

### Input

- `d`  -- direction
- `cp` -- Cartesian product

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
function σ(d::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}
    n1 = dim(cp.X)
    return [σ(d[1:n1], cp.X); σ(d[n1+1:length(d)], cp.Y)]
end

"""
    ρ(d::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}

Return the support function of a Cartesian product.

### Input

- `d`  -- direction
- `cp` -- Cartesian product

### Output

The support function in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
function ρ(d::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}
    n1 = dim(cp.X)
    return ρ(d[1:n1], cp.X) + ρ(d[n1+1:length(d)], cp.Y)
end

"""
    isbounded(cp::CartesianProduct)::Bool

Determine whether a Cartesian product is bounded.

### Input

- `cp` -- Cartesian product

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(cp::CartesianProduct)::Bool
    return isbounded(cp.X) && isbounded(cp.Y)
end

"""
    ∈(x::AbstractVector{N}, cp::CartesianProduct{N})::Bool where {N<:Real}

Check whether a given point is contained in a Cartesian product.

### Input

- `x`  -- point/vector
- `cp` -- Cartesian product

### Output

`true` iff ``x ∈ cp``.
"""
function ∈(x::AbstractVector{N}, cp::CartesianProduct{N})::Bool where {N<:Real}
    @assert length(x) == dim(cp)

    n1 = dim(cp.X)
    return view(x, 1:n1) ∈ cp.X &&
           view(x, n1+1:length(x)) ∈ cp.Y
end

"""
    isempty(cp::CartesianProduct)::Bool

Return if a Cartesian product is empty or not.

### Input

- `cp` -- Cartesian product

### Output

`true` iff any of the sub-blocks is empty.
"""
function isempty(cp::CartesianProduct)::Bool
    return isempty(cp.X) || isempty(cp.Y)
end

"""
    constraints_list(cp::CartesianProduct{N}) where {N<:Real}

Return the list of constraints of a (polyhedral) Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

A list of constraints.
"""
function constraints_list(cp::CartesianProduct{N}) where {N<:Real}
    return constraints_list(CartesianProductArray([cp.X, cp.Y]))
end

"""
    vertices_list(cp::CartesianProduct{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a (polytopic) Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

A list of vertices.

### Algorithm

We assume that the underlying sets are polytopic.
Then the high-dimensional set of vertices is just the Cartesian product of the
low-dimensional sets of vertices.
"""
function vertices_list(cp::CartesianProduct{N}
                      )::Vector{Vector{N}} where {N<:Real}
    # collect low-dimensional vertices lists
    vlist_low = (vertices_list(cp.X), vertices_list(cp.Y))

    # create high-dimensional vertices list
    vlist = Vector{Vector{N}}()
    m = length(vlist_low[1]) * length(vlist_low[2])
    sizehint!(vlist, m)
    for v1 in vlist_low[1]
        for v2 in vlist_low[2]
            push!(vlist, vcat(v1, v2))
        end
    end

    return vlist
end

"""
    linear_map(M::AbstractMatrix{N}, cp::CartesianProduct{N}) where {N<:Real}

Concrete linear map of a (polyhedral) Cartesian product.

### Input

- `M`  -- matrix
- `cp` -- Cartesian product of two convex sets

### Output

A polytope.

### Algorithm

We check if the matrix is invertible.
If so, we convert the Cartesian product to constraint representation.
Otherwise, we convert the Cartesian product to vertex representation.
In both cases, we then call `linear_map` on the resulting polytope.
"""
function linear_map(M::AbstractMatrix{N}, cp::CartesianProduct{N}
                   ) where {N<:Real}
    return linear_map_cartesian_product(M, cp)
end

function linear_map_cartesian_product(M, cp)
    @assert dim(cp) == size(M, 2) "a linear map of size $(size(M)) cannot " *
                                  "be applied to a set of dimension $(dim(cp))"

    if !isinvertible(M)
        # use vertex representation
        P = VPolytope(vertices_list(cp))
    else
        # use constraint representation
        T = isbounded(cp) ? HPolytope : HPolyhedron
        P = T(constraints_list(cp))
    end
    return linear_map(M, P)
end

# ======================================
#  Cartesian product of an array of sets
# ======================================

"""
    CartesianProductArray{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the Cartesian product of a finite number of convex sets.

### Fields

- `array` -- array of sets

### Notes

The `EmptySet` is the absorbing element for `CartesianProductArray`.

Constructors:

- `CartesianProductArray(array::Vector{<:LazySet})` -- default constructor

- `CartesianProductArray([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty product with optional size hint and numeric type
"""
struct CartesianProductArray{N<:Real, S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

# constructor for an empty product with optional size hint and numeric type
function CartesianProductArray(n::Int=0, N::Type=Float64)::CartesianProductArray
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return CartesianProductArray(arr)
end

# EmptySet is the absorbing element for CartesianProductArray
@absorbing(CartesianProductArray, EmptySet)

# add functions connecting CartesianProduct and CartesianProductArray
@declare_array_version(CartesianProduct, CartesianProductArray)

"""
    array(cpa::CartesianProductArray{N, S}
         )::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of a Cartesian product of a finite number of convex sets.

### Input

- `cpa` -- Cartesian product array

### Output

The array of a Cartesian product of a finite number of convex sets.
"""
function array(cpa::CartesianProductArray{N, S}
              )::Vector{S} where {N<:Real, S<:LazySet{N}}
    return cpa.array
end

"""
    dim(cpa::CartesianProductArray)::Int

Return the dimension of a Cartesian product of a finite number of convex sets.

### Input

- `cpa` -- Cartesian product array

### Output

The ambient dimension of the Cartesian product of a finite number of convex
sets.
"""
function dim(cpa::CartesianProductArray)::Int
    return length(cpa.array) == 0 ? 0 : sum([dim(Xi) for Xi in cpa.array])
end

"""
    σ(d::AbstractVector{N}, cpa::CartesianProductArray{N}) where {N<:Real}

Support vector of a Cartesian product array.

### Input

- `d`   -- direction
- `cpa` -- Cartesian product array

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the product sets.
"""
function σ(d::AbstractVector{N}, cpa::CartesianProductArray{N}) where {N<:Real}
    svec = similar(d)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        svec[i0:i1] = σ(d[i0:i1], Xi)
        i0 = i1 + 1
    end
    return svec
end

"""
    ρ(d::AbstractVector{N}, cp::CartesianProductArray{N}) where {N<:Real}

Return the support function of a Cartesian product array.

### Input

- `d`   -- direction
- `cpa` -- Cartesian product array

### Output

The support function in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
function ρ(d::AbstractVector{N}, cpa::CartesianProductArray{N}) where {N<:Real}
    sfun = zero(N)
    i0 = 1
    for Xi in cpa.array
        i1 = i0 + dim(Xi) - 1
        sfun += ρ(d[i0:i1], Xi)
        i0 = i1 + 1
    end
    return sfun
end

"""
    isbounded(cpa::CartesianProductArray)::Bool

Determine whether a Cartesian product of a finite number of convex sets is
bounded.

### Input

- `cpa` -- Cartesian product of a finite number of convex sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cpa::CartesianProductArray)::Bool
    return all(x -> isbounded(x), cpa.array)
end

"""
    ∈(x::AbstractVector{N}, cpa::CartesianProductArray{N}
     )::Bool where {N<:Real}

Check whether a given point is contained in a Cartesian product of a finite
number of sets.

### Input

- `x`   -- point/vector
- `cpa` -- Cartesian product array

### Output

`true` iff ``x ∈ \\text{cpa}``.
"""
function ∈(x::AbstractVector{N}, cpa::CartesianProductArray{N}
          )::Bool where {N<:Real}
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
    isempty(cpa::CartesianProductArray)::Bool

Return if a Cartesian product is empty or not.

### Input

- `cp` -- Cartesian product

### Output

`true` iff any of the sub-blocks is empty.
"""
function isempty(cpa::CartesianProductArray)::Bool
    return any(X -> isempty(X), array(cpa))
end

"""
    constraints_list(cpa::CartesianProductArray{N}) where {N<:Real}

Return the list of constraints of a (polyhedral) Cartesian product of a finite
number of sets.

### Input

- `cpa` -- Cartesian product array

### Output

A list of constraints.
"""
function constraints_list(cpa::CartesianProductArray{N}) where {N<:Real}
    clist = Vector{LinearConstraint{N}}()
    n = dim(cpa)
    sizehint!(clist, n)
    prev_step = 1
    # create high-dimensional constraints list
    for c_low in array(cpa)
        c_low_list = constraints_list(c_low)
        if !isempty(c_low_list)
            indices = prev_step : (dim(c_low_list[1]) + prev_step - 1)
        end
        for constr in c_low_list
            new_constr = LinearConstraint(sparsevec(indices, constr.a, n),
                                          constr.b)
            push!(clist, new_constr)
        end
        prev_step += dim(c_low_list[1])
    end

    return clist
end

"""
    vertices_list(cpa::CartesianProductArray{N}
                 )::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a (polytopic) Cartesian product of a finite
number of sets.

### Input

- `cpa` -- Cartesian product array

### Output

A list of vertices.

### Algorithm

We assume that the underlying sets are polytopic.
Then the high-dimensional set of vertices is just the Cartesian product of the
low-dimensional sets of vertices.
"""
function vertices_list(cpa::CartesianProductArray{N}
                      )::Vector{Vector{N}} where {N<:Real}
    # collect low-dimensional vertices lists
    vlist_low = [vertices_list(X) for X in array(cpa)]

    # create high-dimensional vertices list
    indices_max = [length(vl) for vl in vlist_low]
    m = prod(indices_max)
    vlist = Vector{Vector{N}}(undef, m)
    indices = ones(Int, length(vlist_low))
    v = zeros(N, dim(cpa))
    dim_start_j = 1
    for vl in vlist_low
        v_low = vl[1]
        v[dim_start_j:dim_start_j+length(v_low)-1] = v_low
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
            v[dim_start_j:dim_start_j+length(v_low)-1] = v_low
            dim_start_j += length(v_low)
            j += 1
        end
        indices[j] += 1
        v_low = vlist_low[j][indices[j]]
        v[dim_start_j:dim_start_j+length(v_low)-1] = v_low
    end

    return vlist
end

"""
    same_block_structure(x::AbstractVector{S1}, y::AbstractVector{S2}
                        )::Bool where {S1<:LazySet, S2<:LazySet}

Check whether two vectors of sets have the same block structure, i.e., the
``i``-th entry in the vectors have the same dimension.

### Input

- `x` -- first vector
- `y` -- second vector

### Output

`true` iff the vectors have the same block structure.
"""
function same_block_structure(x::AbstractVector{S1}, y::AbstractVector{S2}
                             )::Bool where {S1<:LazySet, S2<:LazySet}
    if length(x) != length(y)
        return false
    end
    for i in 1:length(x)
        if dim(x[i]) != dim(y[i])
            return false
        end
    end
    return true
end

"""
    block_structure(cpa::CartesianProductArray{N}) where {N}

Returns an array containing the dimension ranges of each block in a `CartesianProductArray`.

### Input

- `cpa` -- Cartesian product array

### Output

A vector of ranges

### Example

```jldoctest
julia> using LazySets: block_structure

julia> cpa = CartesianProductArray([BallInf(zeros(n), 1.0) for n in [3, 1, 2]]);

julia> block_structure(cpa)
3-element Array{UnitRange{Int64},1}:
 1:3
 4:4
 5:6
```
"""
function block_structure(cpa::CartesianProductArray{N}) where {N}
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
    block_to_dimension_indices(cpa::CartesianProductArray{N}, vars::Vector{Int}) where {N}

Returns a vector mapping block index `i`
to tuple `(f, l)` such that either `f = l = -1` or `f` is the first dimension index and `l` is the last dimension
index of the `i`-th block, depending on whether one of the block's dimension indices is specified in `vars`.

### Input

- `cpa` -- Cartesian product array
- `vars` -- list containing the variables of interest, sorted in ascending order

### Output

(i) A vector of tuples, where values in tuple relate to range of dimensions in the i-th block.
(ii) Number of constrained blocks

### Example

```jldoctest
julia> using LazySets: block_to_dimension_indices

julia> cpa = CartesianProductArray([BallInf(zeros(n), 1.0) for n in [1, 3, 2, 3]]);

julia> block_to_dimension_indices(cpa, [2, 4, 8])
(Tuple{Int64,Int64}[(-1, -1), (2, 4), (-1, -1), (7, 9)], 2)
```
This vector represents the mapping "second block from dimension 2 to dimension 4,
fourth block from dimension 7 to dimension 9."
These blocks contain the dimensions specified in `[2, 4, 8]`.
Number of constrained variables here is 2 (2nd and 4th blocks)
"""
function block_to_dimension_indices(cpa::CartesianProductArray{N}, vars::Vector{Int}) where {N}
    result = fill((-1, -1), length(array(cpa)))
    non_empty_length = 0
    start_index, end_index = 1, 0
    v_i = 1
    @inbounds for i in 1:length(cpa.array)
        end_index += dim(cpa.array[i])
        if v_i <= length(vars) && vars[v_i] <= end_index
            result[i] = (start_index, end_index)
            non_empty_length += 1
            while v_i <= length(vars) && vars[v_i] <= end_index
                v_i += 1
            end
        end
        if v_i > length(vars)
            break
        end
        start_index = end_index + 1
    end
    return result, non_empty_length
end

#same but for all variables
function block_to_dimension_indices(cpa::CartesianProductArray{N}) where {N}
    result = Vector{Tuple{Int, Int}}(undef, length(cpa.array))

    start_index, end_index = 1, 0
    @inbounds for i in 1:length(cpa.array)
        end_index += dim(cpa.array[i])
        result[i] = (start_index, end_index)
        start_index = end_index + 1
    end
    return result, length(cpa.array)
end

"""
    substitute_blocks(low_dim_cpa::CartesianProductArray{N},
                         orig_cpa::CartesianProductArray{N},
                           blocks::Vector{Tuple{Int,Int}}) where {N}

Return merged Cartesian Product Array between original CPA and some low-dimensional CPA,
which represents updated subset of variables in specified blocks.

### Input

- `low_dim_cpa` -- low-dimensional cartesian product array
- `orig_cpa` -- original high-dimensional Cartesian product array
- `blocks` -- index of the first variable in each block of `orig_cpa`

### Output

Merged cartesian product array
"""
function substitute_blocks(low_dim_cpa::CartesianProductArray{N},
                           orig_cpa::CartesianProductArray{N},
                           blocks::Vector{Tuple{Int,Int}}) where {N}

    array = Vector{LazySet{N}}(undef, length(orig_cpa.array))
    index = 1
    for bi in 1:length(orig_cpa.array)
        start_ind, end_index = blocks[bi]
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
    linear_map(M::AbstractMatrix{N}, cpa::CartesianProductArray{N}
              ) where {N<:Real}

Concrete linear map of a Cartesian product of a finite number of convex sets.

### Input

- `M`   -- matrix
- `cpa` -- Cartesian product of a finite number of convex sets

### Output

A polytope.

### Algorithm

We check if the matrix is invertible.
If so, we convert the Cartesian product to constraint representation.
Otherwise, we convert the Cartesian product to vertex representation.
In both cases, we then call `linear_map` on the resulting polytope.
"""
function linear_map(M::AbstractMatrix{N}, cpa::CartesianProductArray{N}
                   ) where {N<:Real}
    return linear_map_cartesian_product(M, cpa)
end
