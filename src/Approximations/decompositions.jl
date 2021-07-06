"""
     uniform_partition(n::Int, block_size::Int)

Compute a uniform block partition of the given size.

### Input

- `n`          -- number of dimensions of the partition
- `block_size` -- size of each block

### Output

A vector of ranges, `Vector{UnitRange{Int}}`, such that the size of each block
is the same, if possible.

### Examples

If the number of dimensions `n` is 2, we have two options: either two blocks
of size `1` or one block of size `2`:

```jldoctest partition
julia> LazySets.Approximations.uniform_partition(2, 1)
2-element Vector{UnitRange{Int64}}:
 1:1
 2:2

julia> LazySets.Approximations.uniform_partition(2, 2)
1-element Vector{UnitRange{Int64}}:
 1:2
```

If the block size argument is not compatible with (i.e. does not divide) `n`, the
output is filled with one block of the size needed to reach `n`:

```jldoctest partition
julia> LazySets.Approximations.uniform_partition(3, 1)
3-element Vector{UnitRange{Int64}}:
 1:1
 2:2
 3:3

julia> LazySets.Approximations.uniform_partition(3, 2)
2-element Vector{UnitRange{Int64}}:
 1:2
 3:3

julia> LazySets.Approximations.uniform_partition(10, 6)
2-element Vector{UnitRange{Int64}}:
 1:6
 7:10
```
"""
@inline function uniform_partition(n::Int, block_size::Int)
    m = div(n, block_size)
    r = n % block_size
    res = Vector{UnitRange{Int}}(undef, r > 0 ? m + 1 : m)
    k = 1
    @inbounds for i in 1:m
        l = k + block_size - 1
        res[i] = k:l
        k = l + 1
    end
    if r > 0
        res[m+1] = k:n
    end
    return res
end

"""
    decompose(S::LazySet{N},
              partition::AbstractVector{<:AbstractVector{Int}},
              block_options
             ) where {N}

Decompose a high-dimensional set into a Cartesian product of overapproximations
of the projections over the specified subspaces.

### Input

- `S`             -- set
- `partition`     -- vector of blocks (i.e., of vectors of integers) (see the
                     Notes below)
- `block_options` -- mapping from block indices in `partition` to a
                     corresponding overapproximation option; we only require
                     access via `[⋅]` (but see also the Notes below)

### Output

A `CartesianProductArray` containing the low-dimensional approximated
projections.

### Algorithm

For each block a specific `project` method is called, dispatching on the
corresponding overapproximation option.

### Notes

The argument `partition` requires some discussion.
Typically, the list of blocks should form a partition of the set
``\\{1, \\dots, n\\}`` represented as a list of consecutive blocks, where ``n``
is the ambient dimension of set `S`.

However, technically there is no problem if the blocks are not consecutive,
blocks are missing, blocks occur more than once, or blocks are overlapping.
This function will, however, stick to the order of blocks, so the resulting set
must be interpreted with care in such cases.
One use case is the need of a projection consisting of several blocks.

For convenience, the argument `block_options` can also be given as a single
option instead of a mapping, which is then interpreted as the option for all
blocks.

### Examples

This function supports different options: one can specify the target set, the
degree of accuracy, and template directions.
These options are exemplified below, where we use the following example.

```jldoctest decompose_examples
julia> S = Ball2(zeros(4), 1.0);  # set to be decomposed (4D 2-norm unit ball)

julia> P2d = [1:2, 3:4];  # a partition with two blocks of size two

julia> P1d = [[1], [2], [3], [4]];  # a partition with four blocks of size one
```

#### Different set types

We can decompose using polygons in constraint representation:

```jldoctest decompose_examples
julia> all(ai isa HPolygon for ai in array(decompose(S, P2d, HPolygon)))
true
```

For decomposition into 1D subspaces, we can use `Interval`:

```jldoctest decompose_examples
julia> all(ai isa Interval for ai in array(decompose(S, P1d, Interval)))
true
```

However, if you need to specify different set types for different blocks, the
interface presented so far does not apply.
See the paragraph *Advanced input for different block approximations* below for
how to do that.

#### Refining the decomposition I: ``ε``-close approximation

The ``ε`` option can be used to refine a decomposition, i.e., obtain a more
accurate result.
We use the [Iterative refinement](@ref) algorithm from the `Approximations`
module.

To illustrate this, consider again the set `S` from above.
We decompose into two 2D polygons.
Using smaller ``ε`` implies a better precision, thus more constraints in each 2D
decomposition.
In the following example, we look at the number of constraints in the first
block.

```jldoctest decompose_examples
julia> d(ε, bi) = array(decompose(S, P2d, (HPolygon => ε)))[bi]
d (generic function with 1 method)

julia> [length(constraints_list(d(ε, 1))) for ε in [Inf, 0.1, 0.01]]
3-element Vector{Int64}:
  4
  8
 32
```

#### Refining the decomposition II: template polyhedra

Another way to refine a decomposition is by using template polyhedra.
The idea is to specify a set of template directions and then to compute on each
block the polytopic overapproximation obtained by evaluating the support
function of the given input set over the template directions.

For example, octagonal 2D approximations of the set `S` are obtained with:

```jldoctest decompose_examples
julia> B = decompose(S, P2d, OctDirections);

julia> length(B.array) == 2 && all(dim(bi) == 2 for bi in B.array)
true
```

See [Template directions](@ref) for the available template directions.
Note that, in contrast to the polygonal ``ε``-close approximation from above,
this method can be applied to blocks of any size.

```jldoctest decompose_examples
julia> B = decompose(S, [1:4], OctDirections);

julia> length(B.array) == 1 && dim(B.array[1]) == 4
true
```

#### Advanced input for different block approximations

Instead of defining the approximation option uniformly for each block, we can
define different approximations for different blocks.
The third argument has to be a mapping from block index (in the partition) to
the corresponding approximation option.

For example:

```jldoctest decompose_examples
julia> res = array(decompose(S, P2d, Dict(1 => Hyperrectangle, 2 => 0.1)));

julia> typeof(res[1]), typeof(res[2])
(Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}, HPolygon{Float64, Vector{Float64}})
```
"""
function decompose(S::LazySet{N},
                   partition::AbstractVector{<:AbstractVector{Int}},
                   block_options
                  ) where {N}
    n = dim(S)
    result = Vector{LazySet{N}}(undef, length(partition))

    @inbounds for (i, block) in enumerate(partition)
        result[i] = project(S, block, block_options[i], n)
    end
    return CartesianProductArray(result)
end

# convenience method with uniform block options
function decompose(S::LazySet{N},
                   partition::AbstractVector{<:AbstractVector{Int}},
                   block_options::Union{Type{<:LazySet},
                                        Pair{<:UnionAll, <:Real},
                                        Real,
                                        Type{<:AbstractDirections},
                                        Nothing
                                       }
                  ) where {N}
    n = dim(S)
    result = Vector{LazySet{N}}(undef, length(partition))

    @inbounds for (i, block) in enumerate(partition)
        result[i] = project(S, block, block_options, n)
    end
    return CartesianProductArray(result)
end

# convenience method with uniform block size
function decompose(S::LazySet, block_options; block_size::Int=1)
    partition = uniform_partition(dim(S), block_size)
    return decompose(S, partition, block_options)
end
