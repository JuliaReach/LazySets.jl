export ConvexHull, CH,
       convex_hull!,
       ConvexHull!,
       swap

"""
    ConvexHull{N, S1<:LazySet{N}, S2<:LazySet{N}} <: ConvexSet{N}

Type that represents the convex hull of the union of two sets, i.e., the set

```math
Z = \\{z ∈ ℝ^n : z = λx + (1-λ)y,\\qquad x ∈ X, y ∈ Y,λ ∈ [0, 1] \\}.
```

### Fields

- `X` -- set
- `Y` -- set

### Notes

The `EmptySet` is the neutral element for `ConvexHull`.

This type is always convex.

### Examples

The convex hull of two 100-dimensional Euclidean balls:

```jldoctest
julia> b1, b2 = Ball2(zeros(100), 0.1), Ball2(4*ones(100), 0.2);

julia> c = ConvexHull(b1, b2);

julia> typeof(c)
ConvexHull{Float64, Ball2{Float64, Vector{Float64}}, Ball2{Float64, Vector{Float64}}}
```
"""
struct ConvexHull{N,S1<:LazySet{N},S2<:LazySet{N}} <: ConvexSet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function ConvexHull(X::LazySet{N}, Y::LazySet{N}) where {N}
        @assert dim(X) == dim(Y) "sets in a convex hull must have the same " *
                                 "dimension"
        return new{N,typeof(X),typeof(Y)}(X, Y)
    end
end

isoperationtype(::Type{<:ConvexHull}) = true
concrete_function(::Type{<:ConvexHull}) = convex_hull

isconvextype(::Type{<:ConvexHull}) = true

ispolyhedral(ch::ConvexHull) = ispolyhedral(ch.X) && ispolyhedral(ch.Y)

# EmptySet is the neutral element for ConvexHull
@neutral(ConvexHull, EmptySet)

# Universe is the absorbing element for ConvexHull
@absorbing(ConvexHull, Universe)

# interface for binary set operations
Base.first(ch::ConvexHull) = ch.X
second(ch::ConvexHull) = ch.Y
@declare_binary_operation(ConvexHull)

"""
    CH

Alias for `ConvexHull`.
"""
const CH = ConvexHull

"""
    swap(ch::ConvexHull)

Return a new `ConvexHull` object with the arguments swapped.

### Input

- `ch` -- convex hull of two sets

### Output

A new `ConvexHull` object with the arguments swapped.
"""
function swap(ch::ConvexHull)
    return ConvexHull(ch.Y, ch.X)
end

"""
    dim(ch::ConvexHull)

Return the dimension of a convex hull of two sets.

### Input

- `ch` -- convex hull of two sets

### Output

The ambient dimension of the convex hull of two sets.
"""
function dim(ch::ConvexHull)
    return dim(ch.X)
end

"""
    σ(d::AbstractVector, ch::ConvexHull)

Return a support vector of the convex hull of two sets in a given direction.

### Input

- `d`  -- direction
- `ch` -- convex hull of two sets

### Output

A support vector of the convex hull in the given direction.
"""
function σ(d::AbstractVector, ch::ConvexHull)
    σ1 = σ(d, ch.X)
    σ2 = σ(d, ch.Y)
    ρ1 = dot(d, σ1)
    ρ2 = dot(d, σ2)
    return ρ1 >= ρ2 ? σ1 : σ2
end

"""
    ρ(d::AbstractVector, ch::ConvexHull)

Evaluate the support function of the convex hull of two sets in a given
direction.

### Input

- `d`  -- direction
- `ch` -- convex hull of two sets

### Output

The evaluation of the support function of the convex hull in the given
direction.
"""
function ρ(d::AbstractVector, ch::ConvexHull)
    return max(ρ(d, ch.X), ρ(d, ch.Y))
end

"""
    isbounded(ch::ConvexHull)

Check whether the convex hull of two sets is bounded.

### Input

- `ch` -- convex hull of two sets

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(ch::ConvexHull)
    return isbounded(ch.X) && isbounded(ch.Y)
end

function isboundedtype(::Type{<:ConvexHull{N,S1,S2}}) where {N,S1,S2}
    return isboundedtype(S1) && isboundedtype(S2)
end

"""
    isempty(ch::ConvexHull)

Check whether the convex hull of two sets is empty.

### Input

- `ch` -- convex hull

### Output

`true` iff both wrapped sets are empty.
"""
function isempty(ch::ConvexHull)
    return isempty(ch.X) && isempty(ch.Y)
end

"""
    vertices_list(ch::ConvexHull; [apply_convex_hull]::Bool=true,
                  [backend]=nothing)

Return a list of vertices of the convex hull of two sets.

### Input

- `ch`                -- convex hull of two sets
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)

### Output

A list of vertices.
"""
function vertices_list(ch::ConvexHull;
                       apply_convex_hull::Bool=true,
                       backend=nothing)
    vlist = vcat(vertices_list(ch.X), vertices_list(ch.Y))
    if apply_convex_hull
        convex_hull!(vlist; backend=backend)
    end
    return vlist
end

function translate(ch::ConvexHull, x::AbstractVector)
    X = translate(first(ch), x)
    Y = translate(second(ch), x)
    return ConvexHull(X, Y)
end
