export Zonotope, vertices_list, order

"""
    Zonotope <: LazySet

Type that represents a zonotope.

It is defined as the set

```math
Z = \\left\\{ c + \\sum_{i=1}^p ξ_i g_i,~~ ξ_i \\in [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```
where ``c \\in \\mathbb{R}^n`` is its *center* and ``\\{g_i\\}_{i=1}^p``,
``g_i \\in \\mathbb{R}^n``, is the set of *generators*. This characterization
defines a zonotope as the finite Minkowski sum of line elements. Zonotopes can be
equivalently described as the image of a unit infinity-norm ball in ``\\mathbb{R}^n``
by an affine transformation.

### Fields

- `center`     -- center of the zonotope
- `generators` -- two dimensional matrix where each column is a generator of the zonotope

### Examples

A two-dimensional zonotope with given center and set of generators:

```julia
julia> using LazySets
julia> X = Zonotope([1.0, 0.0], 0.1*eye(2))
LazySets.Zonotope{Float64}([1.0, 0.0], [0.1 0.0; 0.0 0.1])
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
struct Zonotope{N<:Real} <: LazySet
    center::AbstractVector{N}
    generators::AbstractMatrix{N}

    Zonotope{N}(center, generators) where N = new(center, generators)
end
Zonotope(center::AbstractVector{N}, generators::AbstractMatrix{N}) where{N<:Real} = Zonotope{N}(center, generators)


"""
    Zonotope(center::AbstractVector, generators::AbstractVector{T}) where{T<:AbstractVector}

Construct a zonotope given its center and the set of generators.

### Input

- `center`     -- center of the zonotope
- `generators` -- list of generators of the zonotope

### Output

A zonotope with the given center and set of generators.

### Examples

A zonotope in two dimensions with three generators:

```julia
julia> X = Zonotope(ones(2), [[1., 0], [0., 1], [1, 1]])
LazySets.Zonotope{Float64}([1.0, 1.0], [1.0 0.0; 0.0 1.0; 1.0 1.0])
julia> X.generators
2×3 Array{Float64,2}:
 1.0  0.0  1.0
 0.0  1.0  1.0
```
"""
Zonotope(center::AbstractVector, generators::AbstractVector{T}) where{T<:AbstractVector} =
    Zonotope(center, hcat(generators...))

"""
    dim(Z)

Return the ambient dimension of a zonotope.

### Input

- `Z` -- a zonotope

### Output

The ambient dimension of the zonotope.
"""
dim(Z::Zonotope) = length(Z.center)

"""
    vertices_list(Z::Zonotope)

Return the vertices of a zonotope.

### Input

- `Z` -- a zonotope

### Output

The list of vertices as a linear array of ``n``-dimensional vectors, where ``n``
is the dimension of the zonotope.

### Notes

For high dimensions, it would be preferable to develop a `vertex_iterator` approach.
"""
function vertices_list(Z::Zonotope{N})::Vector{Vector{N}} where {N<:Real}
    p = size(Z.generators, 2)
    vlist = Vector{Vector{N}}()
    sizehint!(vlist, 2^p)

    for ξi in IterTools.product([[1, -1] for i = 1:p]...)
        push!(vlist, Z.center .+ Z.generators * collect(ξi))
    end
    # take the convex hull at the output
    return convex_hull!(vlist)
end

"""
    σ(d, Z)

Return the support vector of a zonotope in a given direction.

### Input

- `d` -- direction
- `Z` -- zonotope

### Output

Support vector in the given direction.
"""
function σ(d::AbstractVector, Z::Zonotope)
    return Z.center + Z.generators * unit_step.(Z.generators.' * d)
end

"""
    order(Z::Zonotope)

Return the order of a zonotope.

The order of the zonotope is defined as the quotient between the number of
generators and its dimension.

### Input

- `Z` -- a zonotope

### Output

A rational number representing its order.
"""
order(Z::Zonotope) = size(Z.generators, 2)//dim(Z)
