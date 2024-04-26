import Base: *, ∈, isempty

export CartesianProduct,
       CartesianProduct!,
       swap

"""
    CartesianProduct{N, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the Cartesian product of two sets, i.e., the set

```math
Z = \\{ z ∈ ℝ^{n + m} : z = (x, y),\\qquad x ∈ X, y ∈ Y \\}.
```
If ``X ⊆ ℝ^n`` and ``Y ⊆ ℝ^m``, then ``Z`` is
``n+m``-dimensional.

### Fields

- `X` -- first set
- `Y` -- second set

### Notes

See also [`CartesianProductArray`](@ref) for an implementation of a Cartesian
product of more than two sets.

The `EmptySet` is the almost absorbing element for `CartesianProduct` (except
that the dimension is adapted).

The Cartesian product preserves convexity: if the set arguments are convex, then
their Cartesian product is convex as well.

In some docstrings the word "block" is used to denote each wrapped set, with the
natural order, i.e. we say that the first block of a Cartesian product `cp` is
`cp.X` and the second block is `cp.Y`.

### Examples

The Cartesian product of two sets `X` and `Y` can be constructed either using
`CartesianProduct(X, Y)` or the short-cut notation `X × Y` (to enter the *times*
symbol, write `\\times<tab>`).

```jldoctest cartesianproduct_constructor
julia> I1 = Interval(0, 1);

julia> I2 = Interval(2, 4);

julia> I12 = I1 × I2;

julia> typeof(I12)
CartesianProduct{Float64, Interval{Float64}, Interval{Float64}}
```
A hyperrectangle is the Cartesian product of intervals, so we can convert `I12`
to a `Hyperrectangle` type:

```jldoctest cartesianproduct_constructor
julia> convert(Hyperrectangle, I12)
Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}([0.5, 3.0], [0.5, 1.0])
```
"""
struct CartesianProduct{N,S1<:LazySet{N},S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    function CartesianProduct(X::LazySet{N}, Y::LazySet{N}) where {N}
        return new{N,typeof(X),typeof(Y)}(X, Y)
    end
end

isoperationtype(::Type{<:CartesianProduct}) = true
concrete_function(::Type{<:CartesianProduct}) = cartesian_product

function isconvextype(::Type{CartesianProduct{N,S1,S2}}) where {N,S1,S2}
    return isconvextype(S1) && isconvextype(S2)
end

is_polyhedral(cp::CartesianProduct) = is_polyhedral(cp.X) && is_polyhedral(cp.Y)

# EmptySet is almost the absorbing element for CartesianProduct, but it should
# sum up the dimension of both set arguments
@commutative function CartesianProduct(X::LazySet, ∅::EmptySet)
    N = promote_type(eltype(X), eltype(∅))
    return EmptySet{N}(dim(X) + dim(∅))
end
function CartesianProduct(∅1::EmptySet, ∅2::EmptySet)
    N = promote_type(eltype(∅1), eltype(∅2))
    return EmptySet{N}(dim(∅1) + dim(∅2))
end

# interface for binary set operations
Base.first(cp::CartesianProduct) = cp.X
second(cp::CartesianProduct) = cp.Y
@declare_binary_operation(CartesianProduct)

"""
```
    *(X::LazySet, Y::LazySet)
```

Alias for the binary Cartesian product.
"""
*(X::LazySet, Y::LazySet) = CartesianProduct(X, Y)

"""
    ×(X::LazySet, Y::LazySet)

Alias for the binary Cartesian product.

### Notes

The function symbol can be typed via `\\times<tab>`.
"""
×(X::LazySet, Y::LazySet) = CartesianProduct(X, Y)

"""
    swap(cp::CartesianProduct)

Return a new `CartesianProduct` object with the arguments swapped.

### Input

- `cp` -- Cartesian product

### Output

A new `CartesianProduct` object with the arguments swapped.
"""
function swap(cp::CartesianProduct)
    return CartesianProduct(cp.Y, cp.X)
end

"""
    dim(cp::CartesianProduct)

Return the dimension of a Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

The ambient dimension of the Cartesian product.
"""
function dim(cp::CartesianProduct)
    return dim(cp.X) + dim(cp.Y)
end

"""
    σ(d::AbstractVector, cp::CartesianProduct)

Return a support vector of a Cartesian product.

### Input

- `d`  -- direction
- `cp` -- Cartesian product

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
function σ(d::AbstractVector, cp::CartesianProduct)
    n1 = dim(cp.X)
    return [σ(d[1:n1], cp.X); σ(d[(n1 + 1):length(d)], cp.Y)]
end

# faster version for single-entry vectors
function σ(d::SingleEntryVector, cp::CartesianProduct)
    n1 = dim(cp.X)
    idx = d.i
    if idx <= n1
        return [σ(SingleEntryVector(idx, n1, d.v), cp.X); an_element(cp.Y)]
    else
        n2 = length(d) - n1
        return [an_element(cp.X); σ(SingleEntryVector(idx - n1, n2, d.v), cp.Y)]
    end
end

"""
    ρ(d::AbstractVector, cp::CartesianProduct)

Evaluate the support function of a Cartesian product.

### Input

- `d`  -- direction
- `cp` -- Cartesian product

### Output

The evaluation of the support function in the given direction.
If the direction has norm zero, the result depends on the wrapped sets.
"""
function ρ(d::AbstractVector, cp::CartesianProduct)
    n1 = dim(cp.X)
    return ρ(d[1:n1], cp.X) + ρ(d[(n1 + 1):length(d)], cp.Y)
end

# faster version for single-entry vectors
function ρ(d::SingleEntryVector, cp::CartesianProduct)
    n1 = dim(cp.X)
    idx = d.i
    if idx <= n1
        return ρ(SingleEntryVector(idx, n1, d.v), cp.X)
    else
        n2 = length(d) - n1
        return ρ(SingleEntryVector(idx - n1, n2, d.v), cp.Y)
    end
end

"""
    isbounded(cp::CartesianProduct)

Check whether a Cartesian product is bounded.

### Input

- `cp` -- Cartesian product

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(cp::CartesianProduct)
    return isbounded(cp.X) && isbounded(cp.Y)
end

function isboundedtype(::Type{<:CartesianProduct{N,S1,S2}}) where {N,S1,S2}
    return isboundedtype(S1) && isboundedtype(S2)
end

"""
    ∈(x::AbstractVector, cp::CartesianProduct)

Check whether a given point is contained in a Cartesian product.

### Input

- `x`  -- point/vector
- `cp` -- Cartesian product

### Output

`true` iff ``x ∈ cp``.
"""
function ∈(x::AbstractVector, cp::CartesianProduct)
    @assert length(x) == dim(cp)

    n1 = dim(cp.X)
    return view(x, 1:n1) ∈ cp.X &&
           view(x, (n1 + 1):length(x)) ∈ cp.Y
end

"""
    isempty(cp::CartesianProduct)

Check whether a Cartesian product is empty.

### Input

- `cp` -- Cartesian product

### Output

`true` iff any of the sub-blocks is empty.
"""
function isempty(cp::CartesianProduct)
    return isempty(cp.X) || isempty(cp.Y)
end

"""
    center(cp::CartesianProduct)

Return the center of a Cartesian product of centrally-symmetric sets.

### Input

- `cp` -- Cartesian product of centrally-symmetric sets

### Output

The center of the Cartesian product.
"""
function center(cp::CartesianProduct)
    return vcat(center(cp.X), center(cp.Y))
end

"""
    constraints_list(cp::CartesianProduct)

Return the list of constraints of a (polyhedral) Cartesian product.

### Input

- `cp` -- polyhedral Cartesian product

### Output

A list of constraints.
"""
function constraints_list(cp::CartesianProduct)
    return _constraints_list_cartesian_product(cp)
end

"""
    vertices_list(cp::CartesianProduct)

Return the list of vertices of a (polytopic) Cartesian product.

### Input

- `cp` -- polytopic Cartesian product

### Output

A list of vertices.

### Algorithm

We assume that the underlying sets are polytopic.
Then the high-dimensional set of vertices is just the Cartesian product of the
low-dimensional sets of vertices.
"""
function vertices_list(cp::CartesianProduct)
    # collect low-dimensional vertices lists
    vlist1 = vertices_list(cp.X)
    vlist2 = vertices_list(cp.Y)

    # create high-dimensional vertices list
    N = eltype(cp)
    vlist = Vector{Vector{N}}()
    m = length(vlist1) * length(vlist2)
    sizehint!(vlist, m)
    for v1 in vlist1
        for v2 in vlist2
            push!(vlist, vcat(v1, v2))
        end
    end

    return vlist
end

"""
    linear_map(M::AbstractMatrix, cp::CartesianProduct)

Concrete linear map of a (polyhedral) Cartesian product.

### Input

- `M`  -- matrix
- `cp` -- Cartesian product

### Output

A polytope if `cp` is bounded and a polyhedron otherwise.

### Algorithm

We convert the Cartesian product to constraint representation and then call
`linear_map` on the corresponding polyhedron.

This is a fallback implementation and will fail if the wrapped sets are not
polyhedral.
"""
function linear_map(M::AbstractMatrix, cp::CartesianProduct)
    return _linear_map_cartesian_product(M, cp)
end

function _linear_map_cartesian_product(M, cp)
    @assert dim(cp) == size(M, 2) "a linear map of size $(size(M)) cannot " *
                                  "be applied to a set of dimension $(dim(cp))"

    # use constraint representation
    T = isbounded(cp) ? HPolytope : HPolyhedron
    P = T(constraints_list(cp))
    return linear_map(M, P)
end

function project(cp::CartesianProduct, block::AbstractVector{Int}; kwargs...)
    n1 = dim(cp.X)
    if block[end] <= n1
        # projection completely in the first block
        return project(cp.X, block; kwargs...)
    elseif block[1] > n1
        # projection completely in the second block
        return project(cp.Y, block .- n1; kwargs...)
    end
    # projection is a new Cartesian product of the block-wise projections
    for (i, bi) in enumerate(block)
        if bi > n1
            X = project(cp.X, block[1:(i - 1)]; kwargs...)
            Y = project(cp.Y, block[i:end] .- n1; kwargs...)
            return CartesianProduct(X, Y)
        end
    end
end

"""
    project(cp::CartesianProduct{N, IT, HT}, block::AbstractVector{Int};
            [kwargs...]) where {N, IT<:Interval, HT<:AbstractHyperrectangle{N}}

Concrete projection of the Cartesian product of an interval and a
hyperrectangular set.

### Input

- `cp`    -- Cartesian product of an interval and a hyperrectangle
- `block` -- block structure, a vector with the dimensions of interest

### Output

A hyperrectangle representing the projection of the Cartesian product `cp` on
the dimensions specified by `block`.
"""
function project(cp::CartesianProduct{N,IT,HT}, block::AbstractVector{Int};
                 kwargs...) where {N,IT<:Interval,HT<:AbstractHyperrectangle{N}}
    I = cp.X
    H = cp.Y
    block_vec = collect(block)
    if 1 ∉ block_vec
        block_vec .-= 1
        cH = center(H)
        rH = radius_hyperrectangle(H)
    else
        cH = vcat(center(I), center(H))
        rH = vcat(radius_hyperrectangle(I), radius_hyperrectangle(H))
    end
    return Hyperrectangle(cH[block_vec], rH[block_vec]; check_bounds=false)
end

"""
    project(cp::CartesianProduct{N, IT, ZT}, block::AbstractVector{Int};
            [kwargs...]) where {N, IT<:Interval, ZT<:AbstractZonotope{N}}

Concrete projection of the Cartesian product of an interval and a zonotopic set.

### Input

- `cp`    -- Cartesian product of an interval and a zonotopic set
- `block` -- block structure, a vector with the dimensions of interest

### Output

A zonotope representing the projection of the Cartesian product `cp` on the
dimensions specified by `block`.
"""
function project(cp::CartesianProduct{N,IT,ZT}, block::AbstractVector{Int};
                 kwargs...) where {N,IT<:Interval,ZT<:AbstractZonotope{N}}
    block_vec = collect(block)
    Z = cp.Y
    if 1 ∉ block_vec
        block_vec .-= 1
    else
        Z = convert(Zonotope, cp)
    end
    M = projection_matrix(block_vec, dim(Z), N)
    return linear_map(M, Z)
end

"""
    project(cp::CartesianProduct{N,<:Interval,<:Union{VPolygon,VPolytope}
            block::AbstractVector{Int};
            [kwargs...]) where {N}

Concrete projection of the Cartesian product of an interval and a set in vertex
representation.

### Input

- `cp`    -- Cartesian product of an interval and a `VPolygon` or a `VPolytope`
- `block` -- block structure, a vector with the dimensions of interest

### Output

A `VPolytope` representing the projection of the Cartesian product `cp` on the
dimensions specified by `block`.
"""
function project(cp::CartesianProduct{N,<:Interval,<:Union{VPolygon,VPolytope}},
                 block::AbstractVector{Int};
                 kwargs...) where {N}
    I = cp.X
    P = cp.Y
    block_vec = collect(block)
    if 1 ∉ block_vec
        Pout = project(P, block_vec .- 1; kwargs...)
    else
        out = cartesian_product(I, P)
        Pout = project(out, block_vec; kwargs...)
    end
    return Pout
end

"""
    volume(cp::CartesianProduct)

Compute the volume of a Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

The volume.
"""
function volume(cp::CartesianProduct)
    return volume(cp.X) * volume(cp.Y)
end

function translate(cp::CartesianProduct, x::AbstractVector)
    X = first(cp)
    n1 = dim(X)
    X = translate(X, @view(x[1:n1]))
    Y = translate(second(cp), @view(x[(n1 + 1):end]))
    return CartesianProduct(X, Y)
end
