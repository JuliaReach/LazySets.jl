import Base: ∩

"""
    IntersectionCache

Container for information cached by a lazy `Intersection` object.

### Fields

- `isempty` -- is the intersection empty? There are three possible states,
               encoded as `Int8` values -1, 0, 1:

    * ``-1`` - it is currently unknown whether the intersection is empty
    *  ``0`` - intersection is not empty
    *  ``1`` - intersection is empty
"""
mutable struct IntersectionCache
    isempty::Int8

    # default constructor
    IntersectionCache() = new(Int8(-1))
end

function isempty_known(c::IntersectionCache)
    return c.isempty != Int8(-1)
end

function isempty(c::IntersectionCache)
    @assert isempty_known(c) "'isempty_known' only works if 'isempty' " *
                             "returns 'true'"
    return c.isempty == Int8(1)
end

function set_isempty!(c::IntersectionCache, isempty::Bool)
    return c.isempty = isempty ? Int8(1) : Int8(0)
end

"""
    Intersection{N, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of two sets.

### Fields

- `X`     -- set
- `Y`     -- set
- `cache` -- internal cache for avoiding recomputation; see
             [`IntersectionCache`](@ref)

### Notes

If the arguments of the lazy intersection are half-spaces, the set is simplified
to a polyhedron in constraint representation (`HPolyhedron`).

The intersection preserves convexity: if the set arguments are convex, then
their intersection is convex as well.

The convenience alias `∩` can be typed by `\\cap<tab>`.

### Examples

Create an expression ``Z`` that lazily represents the intersection of two
squares ``X`` and ``Y``:

```jldoctest lazy_intersection
julia> X, Y = BallInf([0.0, 0.0], 0.5), BallInf([1.0, 0.0], 0.75);

julia> Z = X ∩ Y;

julia> typeof(Z)
Intersection{Float64, BallInf{Float64, Vector{Float64}}, BallInf{Float64, Vector{Float64}}}

julia> dim(Z)
2
```

We can check if the intersection is empty with `isempty`:

```jldoctest lazy_intersection
julia> isempty(Z)
false
```

Do not confuse `Intersection` with the concrete operation, which is computed
with the lowercase `intersection` function:

```jldoctest lazy_intersection
julia> W = intersection(X, Y)
Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}([0.375, 0.0], [0.125, 0.5])
```
"""
struct Intersection{N,S1<:LazySet{N},S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2
    cache::IntersectionCache

    # default constructor with dimension check
    function Intersection(X::LazySet{N}, Y::LazySet{N};
                          cache::IntersectionCache=IntersectionCache()) where {N}
        @assert dim(X) == dim(Y) "sets in an intersection must have the same " *
                                 "dimension"
        return new{N,typeof(X),typeof(Y)}(X, Y, cache)
    end
end

# constructors simplifying to HPolyhedron
Intersection(H1::HalfSpace, H2::HalfSpace) = HPolyhedron([H1, H2])
Intersection(H::HalfSpace, P::HPolyhedron) = HPolyhedron(vcat(P.constraints, H))
Intersection(P::HPolyhedron, H::HalfSpace) = HPolyhedron(vcat(P.constraints, H))

∩(X::LazySet, Y::LazySet) = Intersection(X, Y)

concrete_function(::Type{<:Intersection}) = intersection

# Universe is the neutral element for Intersection
@neutral(Intersection, Universe)

# EmptySet is the absorbing element for Intersection
@absorbing(Intersection, EmptySet)

# interface for binary set operations
first(cap::Intersection) = cap.X
second(cap::Intersection) = cap.Y
@declare_binary_operation(Intersection)

"""
    isempty_known(cap::Intersection)

Ask whether the status of emptiness is known.

### Input

- `cap` -- intersection of two sets

### Output

`true` iff the emptiness status is known.
In this case, `isempty(cap)` can be used to obtain the status in constant time.
"""
function isempty_known(cap::Intersection)
    return isempty_known(cap.cache)
end

"""
    set_isempty!(cap::Intersection, isempty::Bool)

Set the status of emptiness in the cache.

### Input

- `cap`     -- intersection of two sets
- `isempty` -- new status of emptiness
"""
function set_isempty!(cap::Intersection, isempty::Bool)
    return set_isempty!(cap.cache, isempty)
end

# equality test ignores the IntersectionCache
function ==(X::Intersection, Y::Intersection)
    if X.X != Y.X
        return false
    end
    if X.Y != Y.Y
        return false
    end
    return true
end

"""
    swap(cap::Intersection)

Return a new `Intersection` object with the arguments swapped.

### Input

- `cap` -- intersection of two sets

### Output

A new `Intersection` object with the arguments swapped.
The old cache is shared between the old and new objects.

### Notes

The advantage of using this function instead of manually swapping the arguments
is that the cache is shared.
"""
function swap(cap::Intersection)
    return Intersection(cap.Y, cap.X; cache=cap.cache)
end

"""
    get_constrained_lowdimset(cpa::CartesianProductArray{N, S},
                              P::AbstractPolyhedron{N}) where {N, S}

Preprocessing step for the intersection between a Cartesian product of a finite
number of sets and a polyhedron.

### Input

- `cpa` -- Cartesian product of a finite number of sets
- `P`   -- polyhedron

### Output

A four-tuple of:
1. a low-dimensional `CartesianProductArray` in the constrained dimensions of
   the original set `cpa`
2. the variables in the constrained blocks,
3. the original block structure of the low-dimensional sets,
4. the list of the constrained blocks.
"""
function get_constrained_lowdimset(cpa::CartesianProductArray{N,S},
                                   P::AbstractPolyhedron{N}) where {N,S}
    if isbounded(P)
        blocks, non_empty_length = block_to_dimension_indices(cpa)
    else
        blocks, non_empty_length = block_to_dimension_indices(cpa, constrained_dimensions(P))
    end

    array = Vector{S}()
    sizehint!(array, non_empty_length)
    variables = Vector{Int}()
    block_structure = Vector{UnitRange{Int}}()
    sizehint!(block_structure, non_empty_length)

    last_var = 1
    for i in eachindex(blocks)
        start_index, end_index = blocks[i]
        block_end = last_var + end_index - start_index
        if start_index != -1
            push!(array, cpa.array[i])
            append!(variables, start_index:end_index)
            push!(block_structure, last_var:block_end)
            last_var = block_end + 1
        end
    end

    return CartesianProductArray(array), variables, block_structure, blocks
end

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
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
