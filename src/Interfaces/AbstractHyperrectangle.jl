import Base: ∈, split
using Base: product

export AbstractHyperrectangle,
       radius_hyperrectangle,
       constraints_list,
       low, high,
       isflat,
       rectify,
       volume

"""
    AbstractHyperrectangle{N<:Real} <: AbstractZonotope{N}

Abstract type for hyperrectangular sets.

### Notes

See [`Hyperrectangle`](@ref) for a standard implementation of this interface.

Every concrete `AbstractHyperrectangle` must define the following functions:
- `radius_hyperrectangle(::AbstractHyperrectangle{N})::Vector{N}` -- return the
    hyperrectangle's radius, which is a full-dimensional vector
- `radius_hyperrectangle(::AbstractHyperrectangle{N}, i::Int)::N` -- return the
    hyperrectangle's radius in the `i`-th dimension
- `isflat(::AbstractHyperrectangle{N})::Bool` -- determine whether the
    hyperrectangle's radius is zero in some dimension

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractHyperrectangle)
5-element Array{Any,1}:
 AbstractSingleton
 BallInf
 Hyperrectangle
 Interval
 SymmetricIntervalHull
```
"""
abstract type AbstractHyperrectangle{N<:Real} <: AbstractZonotope{N} end

isconvextype(::Type{<:AbstractHyperrectangle}) = true

# --- AbstractZonotope interface functions ---


"""
   genmat(H::AbstractHyperrectangle)

Return the generator matrix of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

A matrix where each column represents one generator of `H`.
"""
function genmat(H::AbstractHyperrectangle)
    return genmat_fallback(H)
end

# iterator that wraps the generator matrix
struct HyperrectangleGeneratorIterator{AH<:AbstractHyperrectangle}
    H::AH
    nonflats::Vector{Int}  # dimensions along which `H` is not flat
    dim::Int  # total number of dimensions of `H` (stored for efficiency)

    function HyperrectangleGeneratorIterator(H::AH) where {N<:Real,
            AH<:AbstractHyperrectangle{N}}
        n = dim(H)
        nonflats = Vector{Int}()
        sizehint!(nonflats, n)
        @inbounds for i in 1:n
            if radius_hyperrectangle(H, i) != zero(N)
                push!(nonflats, i)
            end
        end
        return new{AH}(H, nonflats, n)
    end
end

Base.length(it::HyperrectangleGeneratorIterator) = length(it.nonflats)

Base.eltype(::Type{<:HyperrectangleGeneratorIterator{<:AbstractHyperrectangle{N}}}) where {N} =
    SingleEntryVector{N}

function Base.iterate(it::HyperrectangleGeneratorIterator{<:AH},
                      state::Int=1) where {N, AH<:AbstractHyperrectangle{N}}
    if state > length(it.nonflats)
        return nothing
    end
    i = it.nonflats[state]
    r = radius_hyperrectangle(it.H, i)
    g = SingleEntryVector(i, it.dim, r)
    state += 1
    return (g, state)
end

"""
    generators(H::AbstractHyperrectangle)

Return an iterator over the generators of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

An iterator over the generators of `H`.
"""
function generators(H::AbstractHyperrectangle)
    return HyperrectangleGeneratorIterator(H)
end

"""
    ngens(H::AbstractHyperrectangle{N}) where {N<:Real}

Return the number of generators of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

The number of generators.

### Algorithm

A hyperrectangular set has one generator for each non-flat dimension.
"""
function ngens(H::AbstractHyperrectangle{N}) where {N<:Real}
    return sum(i -> radius_hyperrectangle(H, i) > zero(N), 1:dim(H))
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(H::AbstractHyperrectangle{N}
                 )::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

A list of vertices.
Zeros in the radius are correctly handled, i.e., the result does not contain any
duplicate vertices.

### Notes

For high dimensions, it is preferable to develop a `vertex_iterator` approach.

### Algorithm

First we identify the dimensions where `H` is flat, i.e., its radius is zero.
We also compute the number of vertices that we have to create.

Next we create the vertices.
We do this by enumerating all vectors `v` of length `n` (the dimension of `H`)
with entries `-1`/`0`/`1` and construct the corresponding vertex as follows:

```math
    \\text{vertex}(v)(i) = \\begin{cases} c(i) + r(i) & v(i) = 1 \\\\
                                          c(i) & v(i) = 0 \\\\
                                          c(i) - r(i) & v(i) = -1. \\end{cases}
```

For enumerating the vectors `v`, we modify the current `v` from left to right by
changing entries `-1` to `1`, skipping entries `0`, and stopping at the first
entry `1` (but changing it to `-1`).
This way we only need to change the vertex in those dimensions where `v` has
changed, which usually is a smaller number than `n`.
"""
function vertices_list(H::AbstractHyperrectangle{N}
                      )::Vector{Vector{N}} where {N<:Real}
    n = dim(H)

    # identify flat dimensions and store them in a binary vector whose entry in
    # dimension i is 0 if the radius is zero and 1 otherwise
    # the vector will later also contain entries -1
    trivector = Vector{Int8}(undef, n)
    m = 1
    c = center(H)
    v = copy(c)
    @inbounds for i in 1:n
        ri = radius_hyperrectangle(H, i)
        if iszero(ri)
            trivector[i] = Int8(0)
        else
            v[i] += ri
            trivector[i] = Int8(1)
            m *= 2
        end
    end

    # create vertices by modifying the three-valued vector and constructing the
    # corresponding point; for efficiency, we create a copy of the old point and
    # modify every entry that has changed in the three-valued vector
    vlist = Vector{Vector{N}}(undef, m)
    vlist[1] = copy(v)
    @inbounds for i in 2:m
        for j in 1:length(v)
            if trivector[j] == Int8(-1)
                trivector[j] = Int8(1)
                v[j] = c[j] + radius_hyperrectangle(H, j)
            elseif trivector[j] == Int8(1)
                trivector[j] = Int8(-1)
                v[j] = c[j] - radius_hyperrectangle(H, j)
                break
            end
        end
        vlist[i] = copy(v)
    end
    return vlist
end

"""
    constraints_list(H::AbstractHyperrectangle{N}) where {N<:Real}

Return the list of constraints of an axis-aligned hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

A list of linear constraints.
"""
function constraints_list(H::AbstractHyperrectangle{N}) where {N<:Real}
    n = dim(H)
    constraints = Vector{LinearConstraint{N, SingleEntryVector{N}}}(undef, 2*n)
    b, c = high(H), -low(H)
    one_N = one(N)
    for i in 1:n
        ei = SingleEntryVector(i, n, one_N)
        constraints[i] = HalfSpace(ei, b[i])
        constraints[i+n] = HalfSpace(-ei, c[i])
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
    @assert length(d) == dim(H) "a $(length(d))-dimensional vector is " *
                                "incompatible with a $(dim(H))-dimensional set"
    return center(H) .+ sign_cadlag.(d) .* radius_hyperrectangle(H)
end

"""
    ρ(d::AbstractVector{N}, H::AbstractHyperrectangle{N}) where {N<:Real}

Evaluate the support function of a hyperrectangular set in a given direction.

### Input

- `d` -- direction
- `H` -- hyperrectangular set

### Output

Evaluation of the support function in the given direction.
"""
function ρ(d::AbstractVector{N}, H::AbstractHyperrectangle{N}) where {N<:Real}
    @assert length(d) == dim(H) "a $(length(d))-dimensional vector is " *
                                "incompatible with a $(dim(H))-dimensional set"
    c = center(H)
    res = zero(N)
    @inbounds for (i, di) in enumerate(d)
        if di < zero(N)
            res += di * (c[i] - radius_hyperrectangle(H, i))
        elseif di > zero(N)
            res += di * (c[i] + radius_hyperrectangle(H, i))
        end
    end
    return res
end

# helper function for single entry vector
function _ρ_sev_hyperrectangle(d::SingleEntryVector{N}, H::AbstractHyperrectangle{N}) where {N<:Real}

    @assert d.n == dim(H) "a $(d.n)-dimensional vector is " *
                                "incompatible with a $(dim(H))-dimensional set"
    c = center(H)
    return d.v * c[d.i] + abs(d.v) * radius_hyperrectangle(H, d.i)
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
\\max_{x ∈ \\text{vertices}(X)} ‖ x ‖_p = \\max_{α_1, …, α_n ∈ \\{-1, 1\\}} (|c_1 + α_1 r_1|^p + ... + |c_n + α_n r_n|^p)^{1/p}.
```

The function ``x ↦ x^p``, ``p > 0``, is monotonically increasing and thus the
maximum of each term ``|c_i + α_i r_i|^p`` is given by ``|c_i + \\text{sign}(c_i) r_i|^p``
for each ``i``. Hence, ``x^* := \\text{argmax}_{x ∈ X} ‖ x ‖_p`` is the vertex
``c + \\text{diag}(\\text{sign}(c)) r``.
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
    c = center(H)
    @inbounds for i in eachindex(x)
        ri = radius_hyperrectangle(H, i)
        if !_leq(abs(c[i] - x[i]), ri)
            return false
        end
    end
    return true
end

# --- common AbstractHyperrectangle functions ---

"""
    high(H::AbstractHyperrectangle{N})::Vector{N} where {N<:Real}

Return the higher coordinates of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

A vector with the higher coordinates of the hyperrectangular set.
"""
function high(H::AbstractHyperrectangle{N})::Vector{N} where {N<:Real}
    return center(H) .+ radius_hyperrectangle(H)
end

"""
    high(H::AbstractHyperrectangle{N}, i::Int)::N where {N<:Real}

Return the higher coordinate of a hyperrectangular set in a given dimension.

### Input

- `H` -- hyperrectangular set
- `i` -- dimension of interest

### Output

The higher coordinate of the hyperrectangular set in the given dimension.
"""
function high(H::AbstractHyperrectangle{N}, i::Int)::N where {N<:Real}
    return center(H)[i] + radius_hyperrectangle(H, i)
end

"""
    low(H::AbstractHyperrectangle{N})::Vector{N} where {N<:Real}

Return the lower coordinates of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

A vector with the lower coordinates of the hyperrectangular set.
"""
function low(H::AbstractHyperrectangle{N})::Vector{N} where {N<:Real}
    return center(H) .- radius_hyperrectangle(H)
end

"""
    low(H::AbstractHyperrectangle{N}, i::Int)::N where {N<:Real}

Return the lower coordinate of a hyperrectangular set in a given dimension.

### Input

- `H` -- hyperrectangular set
- `i` -- dimension of interest

### Output

The lower coordinate of the hyperrectangular set in the given dimension.
"""
function low(H::AbstractHyperrectangle{N}, i::Int)::N where {N<:Real}
    return center(H)[i] - radius_hyperrectangle(H, i)
end

"""
    isflat(H::AbstractHyperrectangle)::Bool

Determine whether a hyperrectangular set is flat, i.e. whether its radius
is zero in some dimension.

### Input

- `H` -- hyperrectangular set

### Output

`true` iff the hyperrectangular set is flat.

### Notes

For robustness with respect to floating-point inputs, this function relies on
the result of `isapproxzero` when applied to the radius in some dimension.
Hence, this function depends on the absolute zero tolerance `ABSZTOL`.
"""
function isflat(H::AbstractHyperrectangle)::Bool
    return any(i -> isapproxzero(radius_hyperrectangle(H, i)), 1:dim(H))
end

"""
    split(H::AbstractHyperrectangle{N}, num_blocks::AbstractVector{Int}
         ) where {N<:Real}

Partition a hyperrectangular set into uniform sub-hyperrectangles.

### Input

- `H`          -- hyperrectangular set
- `num_blocks` -- number of blocks in the partition for each dimension

### Output

A list of `Hyperrectangle`s.
"""
function split(H::AbstractHyperrectangle{N}, num_blocks::AbstractVector{Int}
              )::Vector{Hyperrectangle{N}} where {N<:Real}
    @assert length(num_blocks) == dim(H) "need number of blocks in each dimension"
    radius = copy(radius_hyperrectangle(H))
    total_number = 1
    lo = low(H)
    hi = high(H)

    # precompute center points in each dimension
    centers = Vector{StepRangeLen{N}}(undef, dim(H))
    for (i, m) in enumerate(num_blocks)
        if m <= 0
            throw(ArgumentError(m, "each dimension needs at least one block"))
        elseif m == one(N)
            centers[i] = range(lo[i] + radius[i], length=1)
        else
            radius[i] /= m
            centers[i] = range(lo[i] + radius[i], step=(2 * radius[i]),
                               length=m)
            total_number *= m
        end
    end

    # create hyperrectangles for every combination of the center points
    result = Vector{Hyperrectangle{N}}(undef, total_number)
    for (i, center) in enumerate(product(centers...))
        result[i] = Hyperrectangle(collect(center), copy(radius))
    end
    return result
end

import LazySets.Arrays: rectify

"""
    rectify(H::AbstractHyperrectangle)

Concrete rectification of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

The `Hyperrectangle` that corresponds to the rectification of `H`.
"""
function rectify(H::AbstractHyperrectangle)
    Hyperrectangle(low=rectify(low(H)), high=rectify(high(H)))
end

"""
    volume(H::AbstractHyperrectangle{N}) where {N<:Real}

Return the volume of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

The volume of ``H``.

### Algorithm

The volume of the ``n``-dimensional hyperrectangle ``H`` with vector radius ``r`` is
``2ⁿ ∏ᵢ rᵢ`` where ``rᵢ`` denotes the ``i``-th component of ``r``.
"""
function volume(H::AbstractHyperrectangle{N}) where {N<:Real}
    vol = mapreduce(x -> 2x, *, radius_hyperrectangle(H))
    return vol
end
