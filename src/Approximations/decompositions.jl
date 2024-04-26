"""
    decompose(S::LazySet{N},
              partition::AbstractVector{<:AbstractVector{Int}},
              block_options) where {N}

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

The argument `partition` deserves some explanation.
Typically, the list of blocks should form a partition of the set
``\\{1, …, n\\}`` represented as a list of consecutive blocks, where ``n``
is the ambient dimension of set `S`.

However, technically there is no problem if the blocks are not consecutive,
blocks are missing, blocks occur more than once, or blocks are overlapping.
The resulting set must be interpreted with care in such cases (e.g., it will not
necessarily be an overapproximation of `S`).

For convenience, the argument `block_options` can also be given as a single
option instead of a mapping, which is then interpreted as the option for all
blocks.

### Examples

The argument `block_options` supports different options: one can specify the
target set, the degree of accuracy, and template directions.
These options are exemplified below, where we use the following example.

```jldoctest decompose_examples
julia> S = Ball2(zeros(4), 1.0);  # set to be decomposed (4D 2-norm unit ball)

julia> P2d = [1:2, 3:4];  # a partition with two blocks, each of size two

julia> P1d = [[1], [2], [3], [4]];  # a partition with four blocks, each of size one
```

#### Different set types

We can decompose using polygons in constraint representation:

```jldoctest decompose_examples
julia> Y = decompose(S, P2d, HPolygon);

julia> all(ai isa HPolygon for ai in Y)
true
```

For decomposition into 1D subspaces, we can use `Interval`:

```jldoctest decompose_examples
julia> Y = decompose(S, P1d, Interval);

julia> all(ai isa Interval for ai in Y)
true
```

However, if you need to specify different set types for different blocks, the
interface presented so far does not apply.
See the paragraph *Advanced input for different block approximations* below for
how to do that.

#### Refining the decomposition: ``ε``-close approximation

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

#### Refining the decomposition: template polyhedra

Another way to refine a decomposition is by using template polyhedra.
The idea is to specify a set of template directions and then compute on each
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
For this purpose, the argument `block_options` can also be a mapping from block
index (in the partition) to the corresponding approximation option.

For example, we can approximate the first block with a `Hyperrectangle` and the
second block with ``ε``-close approximation for ``ε = 0.1``:

```jldoctest decompose_examples
julia> res = array(decompose(S, P2d, Dict(1 => Hyperrectangle, 2 => 0.1)));

julia> typeof(res[1]), typeof(res[2])
(Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}, HPolygon{Float64, Vector{Float64}})
```
"""
function decompose(S::LazySet{N},
                   partition::AbstractVector{<:AbstractVector{Int}},
                   block_options) where {N}
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
                   block_options::Union{Type{<:LazySet},  # set type is not concrete
                                        Pair{<:UnionAll,<:Real},
                                        Real,
                                        Type{<:AbstractDirections},
                                        Nothing}) where {N}
    n = dim(S)
    result = Vector{LazySet{N}}(undef, length(partition))

    @inbounds for (i, block) in enumerate(partition)
        result[i] = project(S, block, block_options, n)
    end
    return CartesianProductArray(result)
end

"""
    decompose(S::LazySet, block_options; [block_size]::Int=1)

Decompose a high-dimensional set into a Cartesian product of overapproximations
of the projections over uniformly-sized subspaces.

### Input

- `S`             -- set
- `block_options` -- overapproximation option or mapping from block indices to a
                     corresponding overapproximation option
- `block_size`    -- (optional; default: `1`) size of the blocks

### Output

A `CartesianProductArray` containing the low-dimensional approximated
projections.
"""
function decompose(S::LazySet, block_options; block_size::Int=1)
    partition = uniform_partition(dim(S), block_size)
    return decompose(S, partition, block_options)
end

# overapproximation of the projection to a fixed target type
function decompose(S::LazySet{N},
                   partition::AbstractVector{<:AbstractVector{Int}},
                   ::Type{ST}) where {N,ST<:LazySet{N}}
    result = Vector{ST}(undef, length(partition))

    @inbounds for (i, block) in enumerate(partition)
        πX = Projection(S, block)
        result[i] = overapproximate(πX, ST)
    end
    return CartesianProductArray(result)
end
