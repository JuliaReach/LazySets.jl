import Base.LinAlg.norm,
       Base.∈

export AbstractHyperrectangle,
       radius_hyperrectangle

"""
    AbstractHyperrectangle{N<:Real} <: AbstractPointSymmetricPolytope{N}

Abstract type for hyperrectangular sets.

### Notes

Every concrete `AbstractHyperrectangle` must define the following functions:
- `radius_hyperrectangle(::AbstractHyperrectangle{N})::Vector{N}` -- return the
    hyperrectangle's radius, which is a full-dimensional vector
- `radius_hyperrectangle(::AbstractHyperrectangle{N}, i::Int)::N` -- return the
    hyperrectangle's radius in the `i`-th dimension

```jldoctest
julia> subtypes(AbstractHyperrectangle)
4-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractSingleton
 LazySets.BallInf
 LazySets.Hyperrectangle
 LazySets.SymmetricIntervalHull
```
"""
abstract type AbstractHyperrectangle{N<:Real} <: AbstractPointSymmetricPolytope{N}
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(H::AbstractHyperrectangle{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

A list of vertices.

### Notes

For high dimensions, it is preferable to develop a `vertex_iterator` approach.
"""
function vertices_list(H::AbstractHyperrectangle{N}
                      )::Vector{Vector{N}} where {N<:Real}
    return [center(H) .+ si .* radius_hyperrectangle(H)
        for si in IterTools.product([[1, -1] for i = 1:dim(H)]...)]
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, H::AbstractHyperrectangle{N}
     )::AbstractVector{N} where {N<:Real}

Return the support vector of a hyperrectangular set in a given direction.

### Input

- `d` -- direction
- `H` -- hyperrectangular set

### Output

The support vector in the given direction.
If the direction has norm zero, the vertex with biggest values is returned.
"""
function σ(d::AbstractVector{N},
           H::AbstractHyperrectangle{N})::AbstractVector{N} where {N<:Real}
    return center(H) .+ sign_cadlag.(d) .* radius_hyperrectangle(H)
end

"""
    norm(H::AbstractHyperrectangle, [p]::Real=Inf)::Real

Return the norm of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.

### Notes

The norm of a hyperrectangular set is defined as the norm of the enclosing ball,
of the given ``p``-norm, of minimal volume.
"""
function norm(H::AbstractHyperrectangle, p::Real=Inf)::Real
    return maximum(map(x -> norm(x, p), vertices_list(H)))
end

"""
    radius(H::AbstractHyperrectangle, [p]::Real=Inf)::Real

Return the radius of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.

### Notes

The radius is defined as the radius of the enclosing ball of the given
``p``-norm of minimal volume with the same center.
It is the same for all corners of a hyperrectangular set.
"""
function radius(H::AbstractHyperrectangle, p::Real=Inf)::Real
    return norm(radius_hyperrectangle(H), p)
end

"""
    diameter(H::AbstractHyperrectangle, [p]::Real=Inf)::Real

Return the diameter of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.

### Notes

The diameter is defined as the maximum distance in the given ``p``-norm between
any two elements of the set.
Equivalently, it is the diameter of the enclosing ball of the given ``p``-norm
of minimal volume with the same center.
"""
function diameter(H::AbstractHyperrectangle, p::Real=Inf)::Real
    return radius(H, p) * 2
end

"""
    ∈(x::AbstractVector{N}, H::AbstractHyperrectangle{N})::Bool where {N<:Real}

Check whether a given point is contained in a hyperrectangular set.

### Input

- `x` -- point/vector
- `H` -- hyperrectangular set

### Output

`true` iff ``x ∈ H``.

### Algorithm

Let ``H`` be an ``n``-dimensional hyperrectangular set, ``c_i`` and ``r_i`` be
the box's center and radius and ``x_i`` be the vector ``x`` in dimension ``i``,
respectively.
Then ``x ∈ H`` iff ``|c_i - x_i| ≤ r_i`` for all ``i=1,…,n``.
"""
function ∈(x::AbstractVector{N},
           H::AbstractHyperrectangle{N})::Bool where {N<:Real}
    @assert length(x) == dim(H)
    for i in eachindex(x)
        if abs(center(H)[i] - x[i]) > radius_hyperrectangle(H, i)
            return false
        end
    end
    return true
end
