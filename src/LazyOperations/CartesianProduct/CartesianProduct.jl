import Base: *

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

concrete_function(::Type{<:CartesianProduct}) = cartesian_product

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
first(cp::CartesianProduct) = cp.X
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

include("center.jl")
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
include("project.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
