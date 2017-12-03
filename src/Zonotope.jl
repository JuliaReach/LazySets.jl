export Zonotope,
       vertices_list,
       order

"""
    Zonotope{N<:Real} <: LazySet

Type that represents a zonotope.

### Fields

- `center`     -- center of the zonotope
- `generators` -- matrix; each column is a generator of the zonotope

### Notes

Mathematically, a zonotope is defined as the set

```math
Z = \\left\\{ c + \\sum_{i=1}^p ξ_i g_i,~~ ξ_i \\in [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```
where ``c \\in \\mathbb{R}^n`` is its *center* and ``\\{g_i\\}_{i=1}^p``,
``g_i \\in \\mathbb{R}^n``, is the set of *generators*. This characterization
defines a zonotope as the finite Minkowski sum of line elements. Zonotopes can
be equivalently described as the image of a unit infinity-norm ball in
``\\mathbb{R}^n`` by an affine transformation.

- `Zonotope(center::AbstractVector{N},
            generators::AbstractMatrix{N}) where {N<:Real}`

- `Zonotope(center::AbstractVector{N},
            generators_list::AbstractVector{T}) where {N<:Real, T<:AbstractVector{N}}`

### Examples

A two-dimensional zonotope with given center and set of generators:

```julia
julia> Z = Zonotope([1.0, 0.0], 0.1*eye(2))
LazySets.Zonotope{Float64}([1.0, 0.0], [0.1 0.0; 0.0 0.1])
julia> dim(Z)
2
```

Compute its vertices:

```julia
julia> vertices_list(Z)
4-element Array{Array{Float64,1},1}:
[0.9, -0.1]
[1.1, -0.1]
[1.1, 0.1]
[0.9, 0.1]
```

Evaluate the support vector in a given direction:

```julia
julia> σ([1., 1.], Z)
2-element Array{Float64,1}:
 1.1
 0.1
```

Alternative constructor: A zonotope in two dimensions with three generators:

```julia
julia> Z = Zonotope(ones(2), [[1., 0.], [0., 1.], [1., 1.]])
LazySets.Zonotope{Float64}([1.0, 1.0], [1.0 0.0; 0.0 1.0; 1.0 1.0])
julia> Z.generators
2×3 Array{Float64,2}:
 1.0  0.0  1.0
 0.0  1.0  1.0
```
"""
struct Zonotope{N<:Real} <: LazySet
    center::AbstractVector{N}
    generators::AbstractMatrix{N}
end
# constructor from center and list of generators
Zonotope(center::AbstractVector{N},
         generators_list::AbstractVector{T}) where {N<:Real, T<:AbstractVector{N}} =
    Zonotope(center, hcat(generators_list...))

"""
    dim(Z::Zonotope)

Return the dimension of a zonotope.

### Input

- `Z` -- zonotope

### Output

The ambient dimension of the zonotope.
"""
dim(Z::Zonotope) = length(Z.center)

"""
    σ(d::AbstractVector{<:Real}, Z::Zonotope)::AbstractVector{<:Real}

Return the support vector of a zonotope in a given direction.

### Input

- `d` -- direction
- `Z` -- zonotope

### Output

Support vector in the given direction.
"""
function σ(d::AbstractVector{<:Real}, Z::Zonotope)::AbstractVector{<:Real}
    return Z.center .+ Z.generators * unit_step.(Z.generators.' * d)
end

"""
    vertices_list(Z::Zonotope{N})::Vector{Vector{N}} where {N<:Real}

Return the vertices of a zonotope.

### Input

- `Z` -- zonotope

### Output

List of vertices.

### Notes

This implementation computes a convex hull.

For high dimensions, it would be preferable to develop a `vertex_iterator`
approach.
"""
function vertices_list(Z::Zonotope{N})::Vector{Vector{N}} where {N<:Real}
    p = size(Z.generators, 2)
    vlist = Vector{Vector{N}}()
    sizehint!(vlist, 2^p)

    for ξi in IterTools.product([[1, -1] for i = 1:p]...)
        push!(vlist, Z.center .+ Z.generators * collect(ξi))
    end

    return convex_hull!(vlist)
end

"""
    order(Z::Zonotope)::Rational

Return the order of a zonotope.

### Input

- `Z` -- zonotope

### Output

A rational number representing the order of the zonotope.

### Notes

The order of a zonotope is defined as the quotient of its number of generators
and its dimension.
"""
order(Z::Zonotope)::Rational = size(Z.generators, 2) // dim(Z)
