"""
    default_block_structure(S::LazySet, set_type::Type{<:LazySet})::AbstractVector{Int}

Compute the default block structure.

### Input

- `S`        -- set
- `set_type` -- target set type

### Output

A vector representing the block structure, such that:

- If the target `set_type` is an interval, the default is blocks of size 1.
- Otherwise, the default is blocks of size 2. Depending on the dimension,
  the last block has size 1 or 2.
"""
@inline function default_block_structure(S::LazySet, set_type::Type{<:LazySet})::AbstractVector{Int}
    n = dim(S)
    if set_type == Interval
        return fill(1, n)
    else
        if n % 2 == 0
            return fill(2, div(n, 2))
        else
            res = fill(2, div(n+1, 2))
            res[end] = 1
            return res
        end
    end
end

"""
    decompose(S::LazySet{N};
              [set_type]::Type{<:Union{HPolygon, Hyperrectangle, Interval}}=Hyperrectangle,
              [ε]::Real=Inf,
              [blocks]::AbstractVector{Int}=default_block_structure(S, set_type),
              [block_types]::Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}(),
              [directions]::Union{Type{<:AbstractDirections}, Void}=nothing
             )::CartesianProductArray where {N<:Real}

Decompose a high-dimensional set into a Cartesian product of overapproximations
of the projections over the specified subspaces.

### Input

- `S`           -- set
- `set_type`    -- (optional, default: `Hyperrectangle`) type of set approximation
                   for each subspace
- `ε`           -- (optional, default: `Inf`) error bound for polytopic approximation
- `blocks`      -- (optional, default: [2, …, 2] or [1, …, 1] if `set_type` is an interval)
                   block structure - a vector with the size of each block
- `block_types` -- (optional, default: Interval for 1D and Hyperrectangle
                   for mD blocks) a mapping from set types to blocks
- `directions`  -- (optional, default: `nothing`) template direction type, or
                   `nothing`

### Output

A `CartesianProductArray` containing the low-dimensional approximated
projections.

### Algorithm

For each block a specific `project` method is called, dispatched on the
`set_type` argument.

### Notes

If `directions` is different from `nothing`, the template directions are used
together with `blocks`.
Otherwise, if `block_types` is given, the options `set_type` and `blocks` are
ignored.

### Examples

The `decompose` function supports different options, such as: supplying different
dimensions for the decomposition, defining the target set of the decomposition,
or specifying the degree of accuracy of the target decomposition. You can also
choose to make the approximations in low dimensions using template directions.
These options are exemplified below.

#### Different dimensions

By default, `decompose` returns a Cartesian product of 2D `Hyperrectangle` sets.
For example:

```jldoctest decompose_examples
julia> import LazySets.Approximations:decompose

julia> S = Ball2(zeros(4), 1.);

julia> array(decompose(S))
2-element Array{LazySets.LazySet{Float64},1}:
 LazySets.Hyperrectangle{Float64}([0.0, 0.0], [1.0, 1.0])
 LazySets.Hyperrectangle{Float64}([0.0, 0.0], [1.0, 1.0])
```

Other block sizes can be specified using the `blocks` option, which refers to
each block size of the partition:

```jldoctest decompose_examples
julia> array(decompose(S, blocks=[1, 3]))
2-element Array{LazySets.LazySet{Float64},1}:
 LazySets.Hyperrectangle{Float64}([0.0], [1.0])
 LazySets.Hyperrectangle{Float64}([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])

julia> array(decompose(S, blocks=[4]))
1-element Array{LazySets.LazySet{Float64},1}:
 LazySets.Hyperrectangle{Float64}([0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0])
```

#### Different set types

We can also decompose using polygons in constraint representation, through the
`set_type` optional argument:

```jldoctest decompose_examples
julia> all([ai isa HPolygon for ai in array(decompose(S, set_type=HPolygon))])
true
```

For decomposition into 1D subspaces, we can use `Interval`:

```jldoctest decompose_examples
julia> all([ai isa Interval for ai in array(decompose(S, set_type=Interval))])
true
```

However, if you need to specify different set types for different blocks, the
interface presented so far does not apply. In the paragraph
*Advanced different set types input* we explain the input `block_types`, that can be
used precisely for that purpose.

#### Refining the decomposition I:  ``ε``-close approximation

The ``ε`` option can be used to refine, that is obtain a more accurate decomposition
in those blocks where `HPolygon` types are used, and it relies on the iterative
refinement algorithm provided in the `Approximations` module.

To illustrate this, consider the unit 4D ball in the 2-norm. Using smaller ``ε``
implies a better precision, thus more constraints in each 2D decomposition:

```jldoctest decompose_examples
julia> S = Ball2(zeros(4), 1.);

julia> d(ε, bi) = array(decompose(S, set_type=HPolygon, ε=ε))[bi]
d (generic function with 1 method)

julia> [length(constraints_list(d(ε, 1))) for ε in [Inf, 0.1, 0.01]]
3-element Array{Int64,1}:
  4
  8
 32
```

#### Refining the decomposition II: template polyhedra

Another way to refine the decomposition is using template polyhedra.
The idea is to specify a set of template directions, and on
each block, compute the polytopic overapproximation obtained by evaluating the
support function of the given input set over the template directions.

For example, octagonal 2D approximations of the ball `S` are obtained with:

```jldoctest decompose_examples
julia> B = decompose(S, directions=OctDirections);

julia> length(B.array) == 2 && all(dim(bi) == 2 for bi in B.array)
true
```

See `template_directions.jl` for the available template directions.
Note that, in contrast to the polygonal ``ε``-close approximation, this method
can be applied for blocks of any size.


```jldoctest decompose_examples
julia> B = decompose(S, directions=OctDirections, blocks=[4]);

julia> length(B.array) == 1 && dim(B.array[1]) == 4
true
```

#### Advanced different set types input

We can define different set types for different blocks, using the
optional `block_types` input argument. It is a dictionary where the keys correspond
to set types, and the values correspond to the blocks, namely the initial and final
block indices should be given.

For example:

```jldoctest decompose_examples
julia> S = Ball2(zeros(3), 1.);

julia> array(decompose(S, block_types=Dict(Interval=>[1:1], Hyperrectangle=>[2:3])))
2-element Array{LazySets.LazySet{Float64},1}:
 LazySets.Interval{Float64,IntervalArithmetic.Interval{Float64}}([-1, 1])
 LazySets.Hyperrectangle{Float64}([0.0, 0.0], [1.0, 1.0])
```

We can additionally pass ε, which is automatically used for each `HPolygon` type block.

```jldoctest decompose_examples
julia> S = Ball2(zeros(8), 1.);

julia> bt = Dict(Interval=>[1:1], Hyperrectangle=>[2:4], HPolygon=>[5:6, 7:8]);

julia> [typeof(ai) for ai in array(decompose(S, block_types=bt, ε=0.01))]
4-element Array{DataType,1}:
 LazySets.Interval{Float64,IntervalArithmetic.Interval{Float64}}
 LazySets.Hyperrectangle{Float64}
 LazySets.HPolygon{Float64}
 LazySets.HPolygon{Float64}
```
"""
function decompose(S::LazySet{N};
                   set_type::Type{<:Union{HPolygon, Hyperrectangle, Interval}}=Hyperrectangle,
                   ε::Real=Inf,
                   blocks::AbstractVector{Int}=default_block_structure(S, set_type),
                   block_types=Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}(),
                   directions::Union{Type{<:AbstractDirections}, Void}=nothing
                  )::CartesianProductArray where {N<:Real}
    n = dim(S)
    result = Vector{LazySet{N}}()

    if directions != nothing
        # template directions
        # potentially defined option set_type is *ignored*       
        block_start = 1
        @inbounds for bi in blocks
            push!(result, project(S, block_start:(block_start + bi - 1), directions(bi), n))
            block_start += bi
        end
    elseif isempty(block_types)
        # use the same target set type for each block
        block_start = 1
        @inbounds for bi in blocks
            push!(result, project(S, block_start:(block_start + bi - 1), set_type, n, ε))
            block_start += bi
        end
    else
        # use potentially different target set type for each block
        # potentially defined options (set_type, blocks) are *ignored*
        initial_block_indices = Vector{Int}()
        blocks = Vector{Int}()
        set_type = Vector{Type{<:LazySet}}()
        @inbounds for (key, val) in block_types
            for bi in val
                push!(set_type, key)
                push!(initial_block_indices, bi[1])
                push!(blocks, bi[end]-bi[1]+1)
            end
        end
        # the second component of this tuple is the starting block index; we
        # assume that blocks do not overlap
        s = sortperm(initial_block_indices)
        blocks = blocks[s]
        set_type = set_type[s]
        block_start = 1
        @inbounds for (i, bi) in enumerate(blocks)
            push!(result, project(S, block_start:(block_start + bi - 1), set_type[i], n, ε))
            block_start += bi
        end
    end
    return CartesianProductArray(result)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:LazySet},
            [n]::Int=dim(S),
            [ε]::Real=Inf
           )::LazySet{N} where {N<:Real}

Default implementation for projecting a high-dimensional set to a given set type
with possible overapproximation.

### Input

- `S` -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type` -- target set type
- `n` -- (optional, default: `dim(S)`) ambient dimension of the set `S`
- `ε` -- (optional, default: `Inf`) ignored

### Output

A set of type `set_type` representing an overapproximation of the projection of
`S`.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block coordinates and zero otherwise.
2. Overapproximate the projected lazy set using `overapproximate`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:LazySet},
                         n::Int=dim(S),
                         ε::Real=Inf
                        )::LazySet{N} where {N<:Real}
    M = sparse(1:length(block), block, ones(N, length(block)), length(block), n)
    return overapproximate(M * S, set_type)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:HPolygon},
            [n]::Int=dim(S),
            [ε]::Real=Inf
           )::HPolygon where {N<:Real}

Project a high-dimensional set to a two-dimensional polygon with a certified
error bound.

### Input

- `S` -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type` -- `HPolygon` - used for dispatch
- `n` -- (optional, default: `dim(S)`) ambient dimension of the set `S`
- `ε` -- (optional, default: `Inf`) error bound for polytopic approximation

### Output

A `HPolygon` representing the epsilon-close approximation of the box
approximation of the projection of `S`.

### Notes

`block` must have length 2.

### Algorithm

If `ε < Inf`, the algorithm proceeds as follows:

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block coordinates and zero otherwise.
2. Overapproximate the set with the given error bound `ε`.

If `ε == Inf`, the algorithm uses a box approximation.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:HPolygon},
                         n::Int=dim(S),
                         ε::Real=Inf
                        )::HPolygon where {N<:Real}
    @assert length(block) == 2 "only 2D HPolygon decomposition is supported"

    # approximation with error bound
    if ε < Inf
        M = sparse([1, 2], [block[1], block[2]], [one(N), one(N)], 2, n)
        return overapproximate(M * S, ε)
    end

    # box approximation
    plus_one = [one(N)]
    minus_one = [-one(N)]
    pe = σ(sparsevec([block[1]], plus_one, n), S)
    pw = σ(sparsevec([block[1]], minus_one, n), S)
    pn = σ(sparsevec([block[2]], plus_one, n), S)
    ps = σ(sparsevec([block[2]], minus_one, n), S)
    pe_bi = dot(DIR_EAST(N), view(pe, block))
    pn_bi = dot(DIR_NORTH(N), view(pn, block))
    pw_bi = dot(DIR_WEST(N), view(pw, block))
    ps_bi = dot(DIR_SOUTH(N), view(ps, block))
    return HPolygon([LinearConstraint(DIR_EAST(N), pe_bi),
                     LinearConstraint(DIR_NORTH(N), pn_bi),
                     LinearConstraint(DIR_WEST(N), pw_bi),
                     LinearConstraint(DIR_SOUTH(N), ps_bi)])
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:Hyperrectangle},
            [n]::Int=dim(S),
            [ε]::Real=Inf
           )::Hyperrectangle where {N<:Real}

Project a high-dimensional set to a low-dimensional hyperrectangle.

### Input

- `S` -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type` -- `Hyperrectangle` - used for dispatch
- `n` -- (optional, default: `dim(S)`) ambient dimension of the set `S`
- `ε` -- (optional, default: `Inf`) - used for dispatch, ignored

### Output

The box approximation of the projection of `S`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:Hyperrectangle},
                         n::Int=dim(S),
                         ε::Real=Inf
                        )::Hyperrectangle where {N<:Real}
    high = Vector{N}(length(block))
    low = similar(high)
    for i in eachindex(block)
        high[i] = σ(sparsevec([block[i]], [one(N)], n), S)[block[i]]
        low[i] = σ(sparsevec([block[i]], [-one(N)], n), S)[block[i]]
    end
    return Hyperrectangle(high=high, low=low)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            directions::AbstractDirections{N},
            n::Int
           )::HPolytope where {N<:Real}

Project a high-dimensional set to a low-dimensional set using template
directions.

### Input

- `S` -- set
- `block` -- block structure - a vector with the dimensions of interest
- `directions` -- template directions
- `n` -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

The template direction approximation of the projection of `S`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         directions::AbstractDirections{N},
                         n::Int=dim(S)
                        )::HPolytope where {N<:Real}
    M = sparse(1:length(block), block, ones(N, length(block)), length(block), n)
    return overapproximate(M * S, directions)
end
