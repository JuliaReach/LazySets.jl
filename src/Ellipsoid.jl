export Ellipsoid

"""
    Ellipsoid <: LazySet

Type that represents an ellipsoid.

It is defined as the set

```math
E = \\left\\{ c + \\sum_{i=1}^p ξ_i g_i,~~ ξ_i \\in [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```
where ``c \\in \\mathbb{R}^n`` is its *center* and ``\\{g_i\\}_{i=1}^p``,
``g_i \\in \\mathbb{R}^n``, is the set of *generators*. This characterization
defines a ellipsoid as the finite Minkowski sum of line elements. ellipsoids can be
equivalently described as the image of a unit infinity-norm ball in ``\\mathbb{R}^n``
by an affine transformation.

### Fields

- `center`     -- center of the ellipsoid
- `generators` -- two dimensional matrix where each column is a generator of the ellipsoid

### Examples

A two-dimensional ellipsoid with given center and set of generators:

```julia
julia> using LazySets
julia> X = ellipsoid([1.0, 0.0], 0.1*eye(2))
LazySets.ellipsoid{Float64}([1.0, 0.0], [0.1 0.0; 0.0 0.1])
julia> dim(X)
2
```

We can ask for its vertices with the `vertices_list` function:

```julia
julia> vertices_list(X)
4-element Array{Array{Float64,1},1}:
[0.9, -0.1]
[1.1, -0.1]
[1.1, 0.1]
[0.9, 0.1]
```

Evaluate the support vector in a given direction:

```julia
julia> σ([1., 1.], X)
2-element Array{Float64,1}:
 1.1
 0.1
```
"""
struct ellipsoid{N<:Real} <: LazySet
    center::AbstractVector{N}
    generators::AbstractMatrix{N}

    ellipsoid{N}(center, generators) where N = new(center, generators)
end
ellipsoid(center::AbstractVector{N}, generators::AbstractMatrix{N}) where{N<:Real} = ellipsoid{N}(center, generators)

"""
    ellipsoid(center::AbstractVector, generators::AbstractVector{T}) where{T<:AbstractVector}

Construct a ellipsoid given its center and the set of generators.

### Input

- `center`     -- center of the ellipsoid
- `generators` -- list of generators of the ellipsoid

### Output

An ellipsoid with the given center and set of generators.

### Examples

A ellipsoid in two dimensions with three generators:

```julia
julia> X = ellipsoid(ones(2), [[1., 0], [0., 1], [1, 1]])
LazySets.ellipsoid{Float64}([1.0, 1.0], [1.0 0.0; 0.0 1.0; 1.0 1.0])
julia> X.generators
2×3 Array{Float64,2}:
 1.0  0.0  1.0
 0.0  1.0  1.0
```
"""
ellipsoid(center::AbstractVector, generators::AbstractVector{T}) where{T<:AbstractVector} =
    ellipsoid(center, hcat(generators...))

"""
    dim(E::ellipsoid)::Int

Return the ambient dimension of an ellipsoid.

### Input

- `E` -- an ellipsoid

### Output

The ambient dimension of the ellipsoid.
"""
dim(E::ellipsoid)::Int = length(E.center)

"""
    σ(d, E)

Return the support vector of an ellipsoid in a given direction.

### Input

- `d` -- direction
- `E` -- ellipsoid

### Output

Support vector in the given direction.
"""
function σ(d::AbstractVector, E::Ellipsoid)
    return Z.center + Z.generators * unit_step.(Z.generators.' * d)
end
