import Base.∈

export AbstractHyperrectangle,
       radius_hyperrectangle,
       constraints_list

"""
    AbstractHyperrectangle{N<:Real} <: AbstractCentrallySymmetricPolytope{N}

Abstract type for hyperrectangular sets.

### Notes

Every concrete `AbstractHyperrectangle` must define the following functions:
- `radius_hyperrectangle(::AbstractHyperrectangle{N})::Vector{N}` -- return the
    hyperrectangle's radius, which is a full-dimensional vector
- `radius_hyperrectangle(::AbstractHyperrectangle{N}, i::Int)::N` -- return the
    hyperrectangle's radius in the `i`-th dimension

```jldoctest
julia> subtypes(AbstractHyperrectangle)
4-element Array{Any,1}:
 AbstractSingleton
 BallInf
 Hyperrectangle
 SymmetricIntervalHull
```
"""
abstract type AbstractHyperrectangle{N<:Real} <: AbstractCentrallySymmetricPolytope{N}
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
    # fast evaluation if H has radius 0
    if radius_hyperrectangle(H) == zeros(N, dim(H))
        return [center(H)]
    end
    return [center(H) .+ si .* radius_hyperrectangle(H)
        for si in Iterators.product([[1, -1] for i = 1:dim(H)]...)][:]
end

"""
    constraints_list(H::AbstractHyperrectangle{N})::Vector{Vector{N}} where {N<:Real}

Return the list of constraints of an axis-aligned hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

A list of linear constraints.
"""
function constraints_list(H::AbstractHyperrectangle{N})::Vector{LinearConstraint{N}} where {N<:Real}
    n = dim(H)
    constraints = Vector{LinearConstraint{N}}(2*n)
    A, b, c = eye(n), high(H), -low(H)
    
    for i in 1:n
        constraints[i] = HalfSpace(A[i, :], b[i])
        constraints[i+n] = HalfSpace(-A[i, :], c[i])
    end
    return constraints
end

# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, H::AbstractHyperrectangle{N}) where {N<:Real}

Return the support vector of a hyperrectangular set in a given direction.

### Input

- `d` -- direction
- `H` -- hyperrectangular set

### Output

The support vector in the given direction.
If the direction has norm zero, the vertex with biggest values is returned.
"""
function σ(d::AbstractVector{N}, H::AbstractHyperrectangle{N}) where {N<:Real}
    return center(H) .+ sign_cadlag.(d) .* radius_hyperrectangle(H)
end

"""
    norm(H::AbstractHyperrectangle, [p]::Real=Inf)::Real

Return the norm of a hyperrectangular set.

The norm of a hyperrectangular set is defined as the norm of the enclosing ball,
of the given ``p``-norm, of minimal volume that is centered in the origin.

### Input

- `H` -- hyperrectangular set
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.

### Algorithm

Recall that the norm is defined as

```math
‖ X ‖ = \\max_{x ∈ X} ‖ x ‖_p = max_{x ∈ \\text{vertices}(X)} ‖ x ‖_p.
```
The last equality holds because the optimum of a convex function over a polytope
is attained at one of its vertices.

This implementation uses the fact that the maximum is achieved in the vertex
``c + \\text{diag}(\\text{sign}(c)) r``, for any ``p``-norm, hence it suffices to
take the ``p``-norm of this particular vertex. This statement is proved below.
Note that, in particular, there is no need to compute the ``p``-norm for *each*
vertex, which can be very expensive. 

If ``X`` is an axis-aligned hyperrectangle and the ``n``-dimensional vectors center
and radius of the hyperrectangle are denoted ``c`` and ``r`` respectively, then
reasoning on the ``2^n`` vertices we have that:

```math
\\max_{x ∈ \\text{vertices}(X)} ‖ x ‖_p = max_{α1, …, αn ∈ {-1, 1}} (|c1 + α1 r1|^p + ... |cn + αn rn|^p)^(1/p).
```

The function ``x ↦ x^p``, ``p > 0``, is monotonically increasing and thus the
maximum of each term ``|ci + αi ri|^p`` is given by ``|ci + sign(ci) ri|^p``
for each ``i``. Hence, ``x^* := argmax_{x ∈ X} ‖ x ‖_p`` is the vertex
``c + diag(sign(c)) r``.
"""
function norm(H::AbstractHyperrectangle, p::Real=Inf)::Real
    c, r = center(H), radius_hyperrectangle(H)
    return norm((@. c + sign_cadlag(c) * r), p)
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
