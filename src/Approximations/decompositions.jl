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
2-element Array{UnitRange{Int64},1}:
 1:1
 2:2

julia> LazySets.Approximations.uniform_partition(2, 2)
1-element Array{UnitRange{Int64},1}:
 1:2
```

If the block size argument is not compatible with (i.e. does not divide) `n`, the
output is filled with one block of the size needed to reach `n`:

```jldoctest partition
julia> LazySets.Approximations.uniform_partition(3, 1)
3-element Array{UnitRange{Int64},1}:
 1:1
 2:2
 3:3

julia> LazySets.Approximations.uniform_partition(3, 2)
2-element Array{UnitRange{Int64},1}:
 1:2
 3:3

julia> LazySets.Approximations.uniform_partition(10, 6)
2-element Array{UnitRange{Int64},1}:
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
             ) where {N<:Real}

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
julia> using LazySets.Approximations: decompose

julia> S = Ball2(zeros(4), 1.);  # set to be decomposed (4D 2-norm unit ball)

julia> P2d = [1:2, 3:4];  # a partition with two blocks of size two

julia> P1d = [[1], [2], [3], [4]];  # a partition with four blocks of size one
```

#### Different set types

We can decompose using polygons in constraint representation:

```jldoctest decompose_examples
julia> all([ai isa HPolygon for ai in array(decompose(S, P2d, HPolygon))])
true
```

For decomposition into 1D subspaces, we can use `Interval`:

```jldoctest decompose_examples
julia> all([ai isa Interval for ai in array(decompose(S, P1d, Interval))])
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
3-element Array{Int64,1}:
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
julia> using LazySets.Approximations: OctDirections

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
(Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}, HPolygon{Float64})
```
"""
function decompose(S::LazySet{N},
                   partition::AbstractVector{<:AbstractVector{Int}},
                   block_options
                  ) where {N<:Real}
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
                  ) where {N<:Real}
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

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            [::Nothing=nothing],
            [n]::Int=dim(S)
           ) where {N<:Real}

Project a high-dimensional set to a given block by using a concrete linear map.

### Input

- `S`       -- set
- `block`   -- block structure - a vector with the dimensions of interest
- `nothing` -- (default: `nothing`) used for dispatch
- `n`       -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set representing the projection of the set `S` to block `block`.

### Algorithm

We apply the function `linear_map`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         ::Nothing=nothing,
                         n::Int=dim(S)
                        ) where {N<:Real}
    m = length(block)
    M = sparse(1:m, block, ones(N, m), m, n)
    return linear_map(M, S)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:LinearMap},
            [n]::Int=dim(S)
           ) where {N<:Real}

Project a high-dimensional set to a given block by using a lazy linear map.

### Input

- `S`         -- set
- `block`     -- block structure - a vector with the dimensions of interest
- `LinearMap` -- used for dispatch
- `n`         -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A lazy `LinearMap` representing the projection of the set `S` to block `block`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:LinearMap},
                         n::Int=dim(S)
                        ) where {N<:Real}
    m = length(block)
    M = sparse(1:m, block, ones(N, m), m, n)
    return M * S
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:LazySet},
            [n]::Int=dim(S)
           ) where {N<:Real}

Project a high-dimensional set to a given block and set type, possibly involving
an overapproximation.

### Input

- `S`        -- set
- `block`    -- block structure - a vector with the dimensions of interest
- `set_type` -- target set type
- `n`        -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set of type `set_type` representing an overapproximation of the projection of
`S`.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected lazy set using `overapproximate` and
`set_type`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:LazySet},
                         n::Int=dim(S)
                        ) where {N<:Real}
    lm = project(S, block, LinearMap, n)
    return overapproximate(lm, set_type)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type_and_precision::Pair{<:UnionAll, <:Real},
            [n]::Int=dim(S)
           ) where {N<:Real}

Project a high-dimensional set to a given block and set type with a certified
error bound.

### Input

- `S`     -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type_and_precision` -- pair `(T, ε)` of a target set type `T` and an
                              error bound `ε` for approximation
- `n`     -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set representing the epsilon-close approximation of the projection of `S`.

### Notes

Currently we only support `HPolygon` as set type, which implies that the set
must be two-dimensional.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected lazy set with the given error bound `ε`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type_and_precision::Pair{<:UnionAll, <:Real},
                         n::Int=dim(S)
                        ) where {N<:Real}
    set_type = set_type_and_precision[1]
    ε = set_type_and_precision[2]
    @assert length(block) == 2 && set_type == HPolygon "currently only 2D " *
        "HPolygon decomposition is supported"

    lm = project(S, block, LinearMap, n)
    return overapproximate(lm, set_type, ε)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            ε::Real,
            [n]::Int=dim(S)
           ) where {N<:Real}

Project a high-dimensional set to a given block and set type with a certified
error bound.

### Input

- `S`     -- set
- `block` -- block structure - a vector with the dimensions of interest
- `ε`     -- error bound for approximation
- `n`     -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set representing the epsilon-close approximation of the projection of `S`.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected lazy set with the given error bound `ε`.
The target set type is chosen automatically.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         ε::Real,
                         n::Int=dim(S)
                        ) where {N<:Real}
    # currently we only support HPolygon
    if length(block) == 2
        set_type = HPolygon
    else
        throw(ArgumentError("ε-close approximation is only supported for 2D " *
                            "blocks"))
    end
    return project(S, block, set_type => ε, n)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            directions::Type{<:AbstractDirections},
            [n]::Int
           ) where {N<:Real}

Project a high-dimensional set to a given block using template directions.

### Input

- `S`          -- set
- `block`      -- block structure - a vector with the dimensions of interest
- `directions` -- template directions
- `n`          -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

The template direction approximation of the projection of `S`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         directions::Type{<:AbstractDirections},
                         n::Int=dim(S)
                        ) where {N<:Real}
    lm = project(S, block, LinearMap, n)
    return overapproximate(lm, directions(length(block)))
end

"""
    project(H::HalfSpace{N}, block::AbstractVector{Int})

Concrete projection of a half-space.

### Input

- `H`        -- set
- `block`    -- block structure, a vector with the dimensions of interest

### Output

A set representing the projection of the half-space `H` on the dimensions
specified by `block`.

### Notes

Currently only the case where the unconstrained dimensions of `H` are a subset
of the `block` variables is implemented.

### Examples

Consider the half-space ``x + y + 0⋅z ≤ 1``, whose ambient dimension is `3`.
The (trivial) projection in the three dimensions is achieved letting the block
of variables to be `[1, 2, 3]`:

```jldoctest project_halfspace
julia> H = HalfSpace([1.0, 1.0, 0.0], 1.0)
HalfSpace{Float64,Array{Float64,1}}([1.0, 1.0, 0.0], 1.0)

julia> using LazySets.Approximations: project

julia> project(H, [1, 2, 3])
HalfSpace{Float64,Array{Float64,1}}([1.0, 1.0, 0.0], 1.0)
```

Projecting along dimensions `1` and `2` only:

```jldoctest project_halfspace
julia> project(H, [1, 2])
HalfSpace{Float64,Array{Float64,1}}([1.0, 1.0], 1.0)
```

In general, use the call syntax `project(H, constrained_dimensions(H))` to return
the half-space projected on the dimensions where it is constrained only:

```jldoctest project_halfspace
julia> project(H, constrained_dimensions(H))
HalfSpace{Float64,Array{Float64,1}}([1.0, 1.0], 1.0)
```
"""
function project(H::HalfSpace{N}, block::AbstractVector{Int}) where {N}
    if constrained_dimensions(H) ⊆ block
        return HalfSpace(H.a[block], H.b)
    else
        error("the concrete projection of a half-space " *
              "for a general block structure is not implemented yet")
    end
end

"""
    project(P::HPolyhedron{N}, block::AbstractVector{Int}) where {N}

Concrete projection of a polyhedron in half-space representation.

### Input

- `P`        -- set
- `block`    -- block structure, a vector with the dimensions of interest

### Output

A set representing the projection of `P` on the dimensions specified by `block`.

### Notes

Currently only the case where the unconstrained dimensions of `P` are a subset
of the `block` variables is implemented.

### Examples

Consider the four-dimensional cross-polytope (unit ball in the 1-norm):

```jldoctest project_hpolyhedron
julia> using LazySets.Approximations: project

julia> P = convert(HPolyhedron, Ball1(zeros(4), 1.0));
```

All dimensions are constrained, and computing the (trivial) projection on the whole
space behaves as expected:

```jldoctest project_hpolyhedron
julia> constrained_dimensions(P)
4-element Array{Int64,1}:
 1
 2
 3
 4

julia> P_1234 = project(P, [1, 2, 3, 4]);

julia> P_1234 == P
true
```
Each constraint of the cross polytope is constrained in all dimensions.

Now let's take a ball in the infinity norm and remove some constraints:

```jldoctest project_hpolyhedron
julia> B = BallInf(zeros(4), 1.0);

julia> c = constraints_list(B)[1:2]
2-element Array{HalfSpace{Float64,LazySets.Arrays.SingleEntryVector{Float64}},1}:
 HalfSpace{Float64,LazySets.Arrays.SingleEntryVector{Float64}}([1.0, 0.0, 0.0, 0.0], 1.0)
 HalfSpace{Float64,LazySets.Arrays.SingleEntryVector{Float64}}([0.0, 1.0, 0.0, 0.0], 1.0)

julia> P = HPolyhedron(c);

julia> constrained_dimensions(P)
2-element Array{Int64,1}:
 1
 2
```

Finally we take the concrete projection onto variables `1` and `2`:

```jldoctest project_hpolyhedron
julia> project(P, [1, 2]) |> constraints_list
2-element Array{HalfSpace{Float64,VN} where VN<:AbstractArray{Float64,1},1}:
 HalfSpace{Float64,Array{Float64,1}}([1.0, 0.0], 1.0)
 HalfSpace{Float64,Array{Float64,1}}([0.0, 1.0], 1.0)
```
"""
function project(P::HPolyhedron{N}, block::AbstractVector{Int}) where {N}
    if constrained_dimensions(P) ⊆ block
        return HPolyhedron([HalfSpace(c.a[block], c.b) for c in constraints_list(P)])
    else
        error("the concrete projection of a polyhedron " *
              "for a general block structure is not implemented yet")
    end
end
