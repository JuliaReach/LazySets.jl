import Base: split
using Base: product

export AbstractHyperrectangle,
       radius_hyperrectangle,
       isflat,
       □

"""
    AbstractHyperrectangle{N} <: AbstractZonotope{N}

Abstract type for hyperrectangular sets.

### Notes

See [`Hyperrectangle`](@ref) for a standard implementation of this interface.

Every concrete `AbstractHyperrectangle` must define the following functions:

- `radius_hyperrectangle(::AbstractHyperrectangle)` -- return the
    hyperrectangle's radius, which is a full-dimensional vector

Among other functions, the following functions are then automatically defined:

- `isflat(::AbstractHyperrectangle)` -- check whether the hyperrectangle's
    radius is zero in some dimension
- `radius_hyperrectangle(::AbstractHyperrectangle, i::Int)` -- return the
    hyperrectangle's radius in the `i`-th dimension

Every hyperrectangular set is also a zonotopic set; see
[`AbstractZonotope`](@ref).

The subtypes of `AbstractHyperrectangle` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractHyperrectangle)
5-element Vector{Any}:
 AbstractSingleton
 BallInf
 Hyperrectangle
 Interval
 SymmetricIntervalHull
```
"""
abstract type AbstractHyperrectangle{N} <: AbstractZonotope{N} end

"""
    radius_hyperrectangle(H::AbstractHyperrectangle)

Return the hyperrectangle radius of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

The hyperrectangle radius of `H`, which is a full-dimensional vector.
"""
function radius_hyperrectangle(::AbstractHyperrectangle) end

"""
    radius_hyperrectangle(H::AbstractHyperrectangle, i::Int)

Return the hyperrectangle radius of a hyperrectangular set in a given dimension.

### Input

- `H` -- hyperrectangular set
- `i` -- dimension

### Output

The hyperrectangle radius of `H` in dimension `i`.

"""
function radius_hyperrectangle(H::AbstractHyperrectangle, i::Int)
    return radius_hyperrectangle(H)[i]
end

"""
    □(c, r)

Convenience constructor of `Hyperrectangle`s or `BallInf`s depending on the type
of `r`.

### Input

- `c` -- center
- `r` -- radius (either a vector for `Hyperrectangle` or a number for `BallInf`)

### Output

A `Hyperrectangle`s or `BallInf`s depending on the type of `r`.

### Notes

The function symbol can be typed via `\\square<tab>`.
"""
function □(c, r) end

function genmat(H::AbstractHyperrectangle)
    gens = generators(H)
    return genmat_fallback(H; gens=gens, ngens=length(gens))
end

# iterator that wraps the generator matrix
struct HyperrectangleGeneratorIterator{AH<:AbstractHyperrectangle}
    H::AH
    nonflats::Vector{Int}  # dimensions along which `H` is not flat
    dim::Int  # total number of dimensions of `H` (stored for efficiency)

    function HyperrectangleGeneratorIterator(H::AH) where {N,
                                                           AH<:AbstractHyperrectangle{N}}
        n = dim(H)
        nonflats = _nonflat_dimensions(H)
        return new{AH}(H, nonflats, n)
    end
end

# return the dimensions of H that are non-flat
function _nonflat_dimensions(H::AbstractHyperrectangle{N}) where {N}
    n = dim(H)
    nonflats = Vector{Int}()
    sizehint!(nonflats, n)
    @inbounds for i in 1:n
        if radius_hyperrectangle(H, i) != zero(N)
            push!(nonflats, i)
        end
    end
    return nonflats
end

Base.length(it::HyperrectangleGeneratorIterator) = length(it.nonflats)

function Base.eltype(::Type{<:HyperrectangleGeneratorIterator{<:AbstractHyperrectangle{N}}}) where {N}
    return SingleEntryVector{N}
end

function Base.iterate(it::HyperrectangleGeneratorIterator{<:AH},
                      state::Int=1) where {N,AH<:AbstractHyperrectangle{N}}
    if state > length(it.nonflats)
        return nothing
    end
    i = it.nonflats[state]
    r = radius_hyperrectangle(it.H, i)
    g = SingleEntryVector(i, it.dim, r)
    state += 1
    return (g, state)
end

function generators(H::AbstractHyperrectangle)
    return HyperrectangleGeneratorIterator(H)
end

"""
# Extended help

    ngens(H::AbstractHyperrectangle)

### Algorithm

A hyperrectangular set has one generator for each non-flat dimension.
"""
function ngens(H::AbstractHyperrectangle)
    N = eltype(H)
    return sum(i -> radius_hyperrectangle(H, i) > zero(N), 1:dim(H))
end

"""
# Extended help

    vertices_list(H::AbstractHyperrectangle; kwargs...)

### Notes

Zeros in the radius are correctly handled, i.e., the result does not contain any
duplicate vertices.

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
function vertices_list(H::AbstractHyperrectangle; kwargs...)
    n = dim(H)

    # identify flat dimensions and store them in a binary vector whose entry in
    # dimension i is 0 if the radius is zero and 1 otherwise
    # the vector will later also contain entries -1
    trivector = Vector{Int8}(undef, n)
    m = 1
    c = center(H)
    v = similar(c)
    copyto!(v, c)
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
    vlist = Vector{typeof(c)}(undef, m)
    vlist[1] = copy(v)
    @inbounds for i in 2:m
        for j in eachindex(v)
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
# Extended help

    constraints_list(H::AbstractHyperrectangle)

### Output

A list of ``2n`` linear constraints, where ``n`` is the dimension of `H`.
"""
function constraints_list(H::AbstractHyperrectangle)
    return _constraints_list_hyperrectangle(H)
end

function _constraints_list_hyperrectangle(H::LazySet{N}) where {N}
    n = dim(H)
    constraints = Vector{HalfSpace{N,SingleEntryVector{N}}}(undef, 2 * n)
    @inbounds for i in 1:n
        ei = SingleEntryVector(i, n, one(N))
        constraints[i] = HalfSpace(ei, high(H, i))
        constraints[i + n] = HalfSpace(-ei, -low(H, i))
    end
    return constraints
end

"""
# Extended help

    σ(d::AbstractVector, H::AbstractHyperrectangle)

### Notes

If the direction vector is zero in dimension ``i``, the result will have the
center's coordinate in that dimension. For instance, for the two-dimensional
infinity-norm ball, if the direction points to the right, the result is the
vector `[1, 0]` in the middle of the right-hand facet.

If the direction has norm zero, the result can be any point in `H`. The default
implementation returns the center of `H`.
"""
@validate function σ(d::AbstractVector, H::AbstractHyperrectangle)
    return center(H) .+ sign_cadlag.(d) .* radius_hyperrectangle(H)
end

# helper function for single-entry vector (used by subtypes)
function _σ_sev_hyperrectangle(d::SingleEntryVector, H::AbstractHyperrectangle)
    N = promote_type(eltype(d), eltype(H))
    s = copy(center(H))
    idx = d.i
    if d.v < zero(N)
        s[idx] -= radius_hyperrectangle(H, idx)
    else
        s[idx] += radius_hyperrectangle(H, idx)
    end
    return s
end

@validate function ρ(d::AbstractVector, H::AbstractHyperrectangle)
    N = promote_type(eltype(d), eltype(H))
    c = center(H)
    res = zero(N)
    @inbounds for (i, di) in enumerate(d)
        ri = radius_hyperrectangle(H, i)
        if di < zero(N)
            res += di * (c[i] - ri)
        elseif di > zero(N)
            res += di * (c[i] + ri)
        end
    end
    return res
end

# helper function for single-entry vector (used by subtypes)
function _ρ_sev_hyperrectangle(d::SingleEntryVector, H::AbstractHyperrectangle)
    return d.v * center(H, d.i) + abs(d.v) * radius_hyperrectangle(H, d.i)
end

"""
# Extended help

    norm(H::AbstractHyperrectangle, [p]::Real=Inf)

### Algorithm

Recall that the norm is defined as

```math
‖ X ‖ = \\max_{x ∈ X} ‖ x ‖_p = max_{x ∈ \\text{vertices}(X)} ‖ x ‖_p.
```
The last equality holds because the optimum of a convex function over a polytope
is attained at one of its vertices.

This implementation uses the fact that the maximum is attained in the vertex
``c + \\text{diag}(\\text{sign}(c)) r`` for any ``p``-norm. Hence it suffices to
take the ``p``-norm of this particular vertex. This statement is proved below.
Note that, in particular, there is no need to compute the ``p``-norm for *each*
vertex, which can be very expensive.

If ``X`` is a hyperrectangle and the ``n``-dimensional vectors center and radius
of the hyperrectangle are denoted ``c`` and ``r`` respectively, then reasoning
on the ``2^n`` vertices we have that:

```math
\\max_{x ∈ \\text{vertices}(X)} ‖ x ‖_p = \\max_{α_1, …, α_n ∈ \\{-1, 1\\}} (|c_1 + α_1 r_1|^p + ... + |c_n + α_n r_n|^p)^{1/p}.
```

The function ``x ↦ x^p``, ``p > 0``, is monotonically increasing and thus the
maximum of each term ``|c_i + α_i r_i|^p`` is given by
``|c_i + \\text{sign}(c_i) r_i|^p`` for each ``i``. Hence,
``x^* := \\text{argmax}_{x ∈ X} ‖ x ‖_p`` is the vertex
``c + \\text{diag}(\\text{sign}(c)) r``.
"""
function norm(H::AbstractHyperrectangle, p::Real=Inf)
    c, r = center(H), radius_hyperrectangle(H)
    return norm((@. c + sign_cadlag(c) * r), p)
end

"""
# Extended help

    radius(H::AbstractHyperrectangle, [p]::Real=Inf)

### Notes

The result is defined as the radius of the enclosing ball of the given
``p``-norm of minimal volume with the same center.
It is the same for all corners of a hyperrectangular set.
"""
function radius(H::AbstractHyperrectangle, p::Real=Inf)
    return norm(radius_hyperrectangle(H), p)
end

"""
# Extended help

    ∈(x::AbstractVector, H::AbstractHyperrectangle)

### Algorithm

Let ``H`` be an ``n``-dimensional hyperrectangular set, ``c_i`` and ``r_i`` be
the center and radius, and ``x_i`` be the vector ``x`` in dimension ``i``,
respectively.
Then ``x ∈ H`` iff ``|c_i - x_i| ≤ r_i`` for all ``i=1,…,n``.
"""
@validate function ∈(x::AbstractVector, H::AbstractHyperrectangle)
    @inbounds for i in eachindex(x)
        ri = radius_hyperrectangle(H, i)
        if !_leq(abs(center(H, i) - x[i]), ri)
            return false
        end
    end
    return true
end

function high(H::AbstractHyperrectangle)
    return center(H) .+ radius_hyperrectangle(H)
end

function high(H::AbstractHyperrectangle, i::Int)
    return center(H, i) + radius_hyperrectangle(H, i)
end

function low(H::AbstractHyperrectangle)
    return center(H) .- radius_hyperrectangle(H)
end

function low(H::AbstractHyperrectangle, i::Int)
    return center(H, i) - radius_hyperrectangle(H, i)
end

function extrema(H::AbstractHyperrectangle)
    c = center(H)
    r = radius_hyperrectangle(H)
    l = c .- r
    h = c .+ r
    return (l, h)
end

function extrema(H::AbstractHyperrectangle, i::Int)
    c = center(H, i)
    r = radius_hyperrectangle(H, i)
    l = c - r
    h = c + r
    return (l, h)
end

"""
    isflat(H::AbstractHyperrectangle)

Check whether a hyperrectangular set is flat, i.e., whether its radius is zero
in some dimension.

### Input

- `H` -- hyperrectangular set

### Output

`true` iff the hyperrectangular set is flat.

### Notes

For robustness with respect to floating-point inputs, this function relies on
the result of `isapproxzero` when applied to the radius in some dimension.
Hence this function depends on the absolute zero tolerance `ABSZTOL`.
"""
function isflat(H::AbstractHyperrectangle)
    return any(i -> isapproxzero(radius_hyperrectangle(H, i)), 1:dim(H))
end

"""
    split(H::AbstractHyperrectangle, num_blocks::AbstractVector{Int})

Partition a hyperrectangular set into uniform sub-hyperrectangles.

### Input

- `H`          -- hyperrectangular set
- `num_blocks` -- number of blocks in the partition for each dimension

### Output

A list of `Hyperrectangle`s.
"""
function split(H::AbstractHyperrectangle, num_blocks::AbstractVector{Int})
    @assert length(num_blocks) == dim(H) "the number of blocks " *
                                         "($(length(num_blocks))) must be specified in each dimension ($(dim(H)))"
    R = radius_hyperrectangle(H)
    T = similar_type(R)
    radius = similar(R)
    copyto!(radius, R)

    total_number = 1
    lo = low(H)
    hi = high(H)
    N = eltype(H)

    # precompute center points in each dimension
    centers = Vector{StepRangeLen{N}}(undef, dim(H))
    @inbounds for (i, m) in enumerate(num_blocks)
        if m <= 0
            throw(ArgumentError("each dimension needs at least one block, got $m"))
        elseif isone(m)
            centers[i] = range(lo[i] + radius[i]; length=1)
        else
            radius[i] /= m
            centers[i] = range(lo[i] + radius[i]; step=(2 * radius[i]),
                               length=m)
            total_number *= m
        end
    end
    radius = convert(T, radius)

    # create hyperrectangles for every combination of the center points
    result = Vector{Hyperrectangle{N,T,T}}(undef, total_number)
    @inbounds for (i, center) in enumerate(product(centers...))
        c = convert(T, collect(center))
        result[i] = Hyperrectangle(c, copy(radius))
    end
    return result
end

# non-uniform partition
function split(H::AbstractHyperrectangle{N},
               bounds::AbstractVector{<:AbstractVector{<:Real}}) where {N}
    @assert length(bounds) == dim(H) "the bounds ($(length(bounds))) must be " *
                                     "specified in each dimension ($(dim(H)))"
    lo = low(H)
    hi = high(H)
    n = length(lo)

    bs = [length(b) + 1 for b in bounds]
    total_number = prod(bs)
    result = Vector{Hyperrectangle{N,Vector{N},Vector{N}}}(undef, total_number)
    for (j, vj) in enumerate(CartesianIterator(bs, false))
        l = Vector{N}(undef, n)
        h = Vector{N}(undef, n)
        for (i, ei) in enumerate(vj)
            if ei == 1
                l[i] = lo[i]
            else
                l[i] = bounds[i][ei - 1]
            end
            if ei == bs[i]
                h[i] = hi[i]
            else
                h[i] = bounds[i][ei]
            end
        end
        result[j] = Hyperrectangle(; low=l, high=h)
    end
    return result
end

"""
# Extended help

    rectify(H::AbstractHyperrectangle)

### Output

A `Hyperrectangle`.
"""
function rectify(H::AbstractHyperrectangle)
    return Hyperrectangle(; low=rectify(low(H)), high=rectify(high(H)))
end

"""
# Extended help

    volume(H::AbstractHyperrectangle)

### Algorithm

The volume of the ``n``-dimensional hyperrectangle ``H`` with radius vector
``r`` is ``2ⁿ ∏ᵢ rᵢ`` where ``rᵢ`` denotes the ``i``-th component of ``r``.
"""
function volume(H::AbstractHyperrectangle)
    return _volume_hyperrectangle(H)
end

function _volume_hyperrectangle(H)
    return mapreduce(x -> 2x, *, radius_hyperrectangle(H))
end

function _surface_area_hyperrectangle(H)
    l, b, h = radius_hyperrectangle(H)
    l *= 2
    b *= 2
    h *= 2
    return 2 * (l * b + l * h + b * h)
end

@validate function area(H::AbstractHyperrectangle)
    if dim(H) == 2
        return _volume_hyperrectangle(H)
    else
        return _surface_area_hyperrectangle(H)
    end
end

@validate function project(H::AbstractHyperrectangle, block::AbstractVector{Int}; kwargs...)
    πc = center(H)[block]
    πr = radius_hyperrectangle(H)[block]
    return Hyperrectangle(πc, πr; check_bounds=false)
end

@validate_commutative function distance(x::AbstractVector, H::AbstractHyperrectangle; p::Real=2)
    # compute closest point
    y = similar(x)
    outside = false
    @inbounds for i in eachindex(x)
        ci = center(H, i)
        ri = radius_hyperrectangle(H, i)
        d = x[i] - ci
        if abs(d) <= ri
            # point is inside in the projection → y[i] is x[i]
            y[i] = x[i]
        else
            # point is outside in the projection → y[i] is on the border
            y[i] = ci + sign_cadlag(d) * ri
            outside = true
        end
    end

    if !outside
        # point is inside
        N = promote_type(eltype(x), eltype(H))
        return zero(N)
    end

    return distance(x, y; p=p)
end

"""
# Extended help

    reflect(H::AbstractHyperrectangle)

### Output

A `Hyperrectangle`.

### Algorithm

If ``H`` has center ``c`` and radius ``r``, then ``-H`` has center ``-c`` and
radius ``r``.
"""
function reflect(H::AbstractHyperrectangle)
    return Hyperrectangle(-center(H), radius_hyperrectangle(H))
end
