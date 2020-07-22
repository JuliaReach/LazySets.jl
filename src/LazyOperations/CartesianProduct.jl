import Base: *, ∈, isempty

export CartesianProduct,
       CartesianProduct!,
       swap

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
Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}([0.5, 3.0], [0.5, 1.0])
```
"""
struct CartesianProduct{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2
end

isoperationtype(::Type{<:CartesianProduct}) = true
isconvextype(::Type{CartesianProduct{N, S1, S2}}) where {N, S1, S2} = isconvextype(S1) && isconvextype(S2)

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
    isbounded(cp::CartesianProduct)

Determine whether a Cartesian product is bounded.

### Input

- `cp` -- Cartesian product

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(cp::CartesianProduct)
    return isbounded(cp.X) && isbounded(cp.Y)
end

"""
    ∈(x::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}

Check whether a given point is contained in a Cartesian product.

### Input

- `x`  -- point/vector
- `cp` -- Cartesian product

### Output

`true` iff ``x ∈ cp``.
"""
function ∈(x::AbstractVector{N}, cp::CartesianProduct{N}) where {N<:Real}
    @assert length(x) == dim(cp)

    n1 = dim(cp.X)
    return view(x, 1:n1) ∈ cp.X &&
           view(x, n1+1:length(x)) ∈ cp.Y
end

"""
    isempty(cp::CartesianProduct)

Return if a Cartesian product is empty or not.

### Input

- `cp` -- Cartesian product

### Output

`true` iff any of the sub-blocks is empty.
"""
function isempty(cp::CartesianProduct)
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
    vertices_list(cp::CartesianProduct{N}) where {N<:Real}

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
function vertices_list(cp::CartesianProduct{N}) where {N<:Real}
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

A polytope if `cp` is bounded and a polyhedron otherwise.

### Algorithm

We convert the Cartesian product to constraint representation and then call
`linear_map` on the corresponding polyhedron.
"""
function linear_map(M::AbstractMatrix{N}, cp::CartesianProduct{N}
                   ) where {N<:Real}
    return linear_map_cartesian_product(M, cp)
end

function linear_map_cartesian_product(M, cp)
    @assert dim(cp) == size(M, 2) "a linear map of size $(size(M)) cannot " *
                                  "be applied to a set of dimension $(dim(cp))"

    # use constraint representation
    T = isbounded(cp) ? HPolytope : HPolyhedron
    P = T(constraints_list(cp))
    return linear_map(M, P)
end

function concretize(cp::CartesianProduct)
    return cartesian_product(concretize(cp.X), concretize(cp.Y))
end
