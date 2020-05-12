import Base: *, ∈, isempty

export CartesianProductArray,
       array,
       same_block_structure

# =======================================
#  Cartesian product of an array of sets
# =======================================

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
function CartesianProductArray(n::Int=0, N::Type=Float64)
   arr = Vector{LazySet{N}}()
   sizehint!(arr, n)
   return CartesianProductArray(arr)
end

isoperationtype(::Type{<:CartesianProductArray}) = true
isconvextype(::Type{CartesianProductArray{N, S}}) where {N, S} = isconvextype(S)

# EmptySet is the absorbing element for CartesianProductArray
@absorbing(CartesianProductArray, EmptySet)

# add functions connecting CartesianProduct and CartesianProductArray
@declare_array_version(CartesianProduct, CartesianProductArray)

"""
   array(cpa::CartesianProductArray{N, S}) where {N<:Real, S<:LazySet{N}}

Return the array of a Cartesian product of a finite number of convex sets.

### Input

- `cpa` -- Cartesian product array

### Output

The array of a Cartesian product of a finite number of convex sets.
"""
function array(cpa::CartesianProductArray{N, S}
             ) where {N<:Real, S<:LazySet{N}}
   return cpa.array
end

"""
   dim(cpa::CartesianProductArray)

Return the dimension of a Cartesian product of a finite number of convex sets.

### Input

- `cpa` -- Cartesian product array

### Output

The ambient dimension of the Cartesian product of a finite number of convex
sets.
"""
function dim(cpa::CartesianProductArray)
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

# faster version for sparse vectors
function σ(d::AbstractSparseVector{N}, cpa::CartesianProductArray{N}
         ) where {N<:Real}
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

# faster version for sparse vectors
function ρ(d::AbstractSparseVector{N}, cpa::CartesianProductArray{N}
         ) where {N<:Real}
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

"""
   isbounded(cpa::CartesianProductArray)

Determine whether a Cartesian product of a finite number of convex sets is
bounded.

### Input

- `cpa` -- Cartesian product of a finite number of convex sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cpa::CartesianProductArray)
   return all(isbounded, cpa.array)
end

"""
   ∈(x::AbstractVector{N}, cpa::CartesianProductArray{N}) where {N<:Real}

Check whether a given point is contained in a Cartesian product of a finite
number of sets.

### Input

- `x`   -- point/vector
- `cpa` -- Cartesian product array

### Output

`true` iff ``x ∈ \\text{cpa}``.
"""
function ∈(x::AbstractVector{N}, cpa::CartesianProductArray{N}
         ) where {N<:Real}
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

Return if a Cartesian product is empty or not.

### Input

- `cp` -- Cartesian product

### Output

`true` iff any of the sub-blocks is empty.
"""
function isempty(cpa::CartesianProductArray)
   return any(isempty, array(cpa))
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
   clist = Vector{LinearConstraint{N, SparseVector{N, Int}}}()
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
           new_constr = LinearConstraint(sparsevec(indices, constr.a, n), constr.b)
           push!(clist, new_constr)
       end
       prev_step += dim(c_low_list[1])
   end

   return clist
end

"""
   vertices_list(cpa::CartesianProductArray{N}) where {N<:Real}

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
function vertices_list(cpa::CartesianProductArray{N}) where {N<:Real}
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
                       ) where {S1<:LazySet, S2<:LazySet}

Check whether two vectors of sets have the same block structure, i.e., the
``i``-th entry in the vectors have the same dimension.

### Input

- `x` -- first vector
- `y` -- second vector

### Output

`true` iff the vectors have the same block structure.
"""
function same_block_structure(x::AbstractVector{S1}, y::AbstractVector{S2}
                            ) where {S1<:LazySet, S2<:LazySet}
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

julia> m, k = block_to_dimension_indices(cpa, [2, 4, 8]);

julia> m
4-element Array{Tuple{Int64,Int64},1}:
 (-1, -1)
 (2, 4)
 (-1, -1)
 (7, 9)

julia> k
2
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
