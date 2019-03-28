import LazySets.dim
using LazySets: isapproxzero

"""
    AbstractDirections{N}

Abstract type for template direction representations.

### Notes

All subtypes should implement the standard iterator methods from `Base` and the
function `dim(d<:AbstractDirections)::Int`.
"""
abstract type AbstractDirections{N} end

"""
    dim(ad::AbstractDirections)::Int

Returns the dimension of the generated directions.

### Input

- `ad` -- box direction representation

### Output

The dimension of the generated directions.
"""
function dim(ad::AbstractDirections)::Int
    return ad.n
end

# ==================================================
# Box directions
# ==================================================

"""
    UnitVector{T} <: AbstractVector{T}

A lazy unit vector with arbitrary one-element.

### Fields

- `i` -- index of non-zero entry
- `n` -- vector length
- `v` -- non-zero entry
"""
struct UnitVector{T} <: AbstractVector{T}
    i::Int
    n::Int
    v::T
end

function Base.getindex(e::UnitVector{T}, i::Int) where T
    @boundscheck @assert 1 <= i <= e.n
    return i == e.i ? e.v : zero(T)
end

Base.size(e::UnitVector) = (e.n,)

"""
    BoxDirections{N} <: AbstractDirections{N}

Box direction representation.

### Fields

- `n` -- dimension
"""
struct BoxDirections{N} <: AbstractDirections{N}
    n::Int
end

# constructor for type Float64
BoxDirections(n::Int) = BoxDirections{Float64}(n)

Base.eltype(::Type{BoxDirections{N}}) where {N} = AbstractVector{N}
Base.length(bd::BoxDirections) = 2 * bd.n

@static if VERSION < v"0.7-"
@eval begin

Base.start(bd::BoxDirections) = 1
Base.next(bd::BoxDirections{N}, state) where {N} = (
    UnitVector{N}(abs(state), bd.n, convert(N, sign(state))), # value
    state == bd.n ? -bd.n : state + 1) # next state
Base.done(bd::BoxDirections, state) = state == 0

end # @eval
else
@eval begin

function Base.iterate(bd::BoxDirections{N}, state::Int=1) where {N}
    if state == 0
        return nothing
    end
    vec = UnitVector{N}(abs(state), bd.n, convert(N, sign(state)))
    state = (state == bd.n) ? -bd.n : state + 1
    return (vec, state)
end

end # @eval
end # if

# ==================================================
# Octagonal directions
# ==================================================

"""
    OctDirections{N} <: AbstractDirections{N}

Octagon direction representation.

### Fields

- `n` -- dimension

### Notes

Octagon directions consist of all vectors that are zero almost everywhere except
in two dimensions ``i``, ``j`` (possibly ``i = j``) where it is ``±1``.
"""
struct OctDirections{N} <: AbstractDirections{N}
    n::Int
end

# constructor for type Float64
OctDirections(n::Int) = OctDirections{Float64}(n)

Base.eltype(::Type{OctDirections{N}}) where {N} = AbstractVector{N}
Base.length(od::OctDirections) = 2 * od.n^2

@static if VERSION < v"0.7-"
@eval begin

function Base.start(od::OctDirections{N}) where {N}
    if od.n == 1
        return 1 # fall back to box directions in 1D case
    end
    vec = zeros(N, od.n)
    vec[1] = one(N)
    vec[2] = vec[1]
    return (vec, 1, 2)
end

function Base.next(od::OctDirections{N}, state::Tuple) where {N}
    vec = state[1]
    i = state[2]
    j = state[3]
    if vec[j] > 0
        # change sign at j
        vec[j] = -one(N)
    elseif vec[i] > 0
        # change sign at i and j
        vec[i] = -one(N)
        vec[j] = one(N)
    elseif j < od.n
        # advance j, reset i
        vec[i] = one(N)
        vec[j] = zero(N)
        j = j + 1
        vec[j] = one(N)
    elseif i < od.n - 1
        # advance i
        vec[i] = zero(N)
        vec[j] = vec[i]
        i = i + 1
        j = i + 1
        vec[i] = one(N)
        vec[j] = vec[i]
    else
        vec = zeros(N, od.n)
        vec[1] = one(N)
        vec[2] = vec[1]
        return (vec, 1) # continue with box directions
    end
    return (copy(vec), (vec, i, j))
end

Base.next(od::OctDirections{N}, state::Int) where {N} =
    next(BoxDirections{N}(od.n), state)

Base.done(od::OctDirections, state) = state == 0

end # @eval
else
@eval begin

function Base.iterate(od::OctDirections{N}) where {N}
    if od.n == 1
        # fall back to box directions in 1D case
        return iterate(od, 1)
    end
    vec = spzeros(N, od.n)
    vec[1] = one(N)
    vec[2] = one(N)
    return (copy(vec), (vec, 1, 2))
end

# basic idea: modify the vector from the previous iteration
# have two indices i = 1, j = i + 1
# - run j through all indices > i
# - i = i + 1, j = i + 1
# - for any pair (i, j) create four vectors [..., i: ±1, ..., j: ±1, ...]
# in the end continue with box directions
function Base.iterate(od::OctDirections{N}, state::Tuple) where {N}
    # continue with octagon directions
    vec = state[1]
    i = state[2]
    j = state[3]
    if vec[j] > 0
        # change sign at j
        vec[j] = -one(N)
    elseif vec[i] > 0
        # change sign at i and j
        vec[i] = -one(N)
        vec[j] = one(N)
    elseif j < od.n
        # advance j, reset i
        vec[i] = one(N)
        vec[j] = zero(N)
        j = j + 1
        vec[j] = one(N)
    elseif i < od.n - 1
        # advance i
        vec[i] = zero(N)
        vec[j] = zero(N)
        i = i + 1
        j = i + 1
        vec[i] = one(N)
        vec[j] = one(N)
    else
        # continue with box directions
        return iterate(od, 1)
    end
    return (copy(vec), (vec, i, j))
end

function Base.iterate(od::OctDirections{N}, state::Int) where {N}
    # continue with box directions
    return iterate(BoxDirections{N}(od.n), state)
end

end # @eval
end # if

# ==================================================
# Box-diagonal directions
# ==================================================

"""
    BoxDiagDirections{N} <: AbstractDirections{N}

Box-diagonal direction representation.

### Fields

- `n` -- dimension

### Notes

Box-diagonal directions can be seen as the union of diagonal directions (all
entries are ±1) and box directions (one entry is ±1, all other entries are 0).
The iterator first enumerates all diagonal directions, and then all box
directions.
"""
struct BoxDiagDirections{N} <: AbstractDirections{N}
    n::Int
end

# constructor for type Float64
BoxDiagDirections(n::Int) = BoxDiagDirections{Float64}(n)

Base.eltype(::Type{BoxDiagDirections{N}}) where {N} = AbstractVector{N}
Base.length(bdd::BoxDiagDirections) = bdd.n == 1 ? 2 : 2^bdd.n + 2 * bdd.n

@static if VERSION < v"0.7-"
@eval begin

Base.start(bdd::BoxDiagDirections{N}) where {N} = ones(N, bdd.n)

function Base.next(bdd::BoxDiagDirections{N}, state::AbstractVector) where {N}
    i = 1
    while i <= bdd.n && state[i] < 0
        state[i] = -state[i]
        i = i+1
    end
    if i > bdd.n
        if bdd.n == 1
            return (copy(state), 0) # finish here to avoid duplicates
        else
            return (copy(state), 1) # continue with box directions
        end
    else
        state[i] = -state[i]
        return (copy(state), state)
    end
end

Base.next(bdd::BoxDiagDirections{N}, state::Int) where {N} =
    next(BoxDirections{N}(bdd.n), state)
Base.done(bdd::BoxDiagDirections, state) = state == 0

end # @eval
else
@eval begin

Base.iterate(bdd::BoxDiagDirections{N}) where {N} =
    (ones(N, bdd.n), ones(N, bdd.n))

function Base.iterate(bdd::BoxDiagDirections{N},
                      state::AbstractVector) where {N}
    # continue with diagonal directions
    i = 1
    while i <= bdd.n && state[i] < 0
        state[i] = -state[i]
        i = i+1
    end
    if i > bdd.n
        if bdd.n == 1
            # finish here to avoid duplicates
            return nothing
        else
            # continue with box directions
            return iterate(bdd, 1)
        end
    else
        state[i] = -state[i]
        return (copy(state), state)
    end
end

function Base.iterate(bdd::BoxDiagDirections{N}, state::Int) where {N}
    # continue with box directions
    return iterate(BoxDirections{N}(bdd.n), state)
end

end # @eval
end # if

# ==================================================
# Spherical directions
# ==================================================

# returns `true` if and only if there is a vector in the list of vectors S
# such that it is approximately d
function _isapprox_included(d::Vector{N}, S::Vector{Vector{N}}) where {N}
    result = false
    for si in S
        if isapproxzero(norm(d-si))
            result = true
            break
        end
    end
    return result
end

"""
    SphericalDirections{N} <: AbstractDirections{N}

Spherical directions representation.

### Fields

- `Nθ`    -- length of the partition of the azimuthal angle
- `Nφ`    -- length of the partition of the polar angle
- `stack` -- list of computed directions

### Notes

The `SphericalDirections` constructor provides a sample of the unit sphere
in ``\\mathbb{R}^3``, which is parameterized by the azimuthal and polar angles
``θ ∈ Dθ := [-π, π]`` and ``φ ∈ Dφ := [0, 2π]`` respectively, see the wikipedia
entry [Spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system).
The domains ``Dθ`` and ``Dφ`` are discretized in ``Nθ`` and ``Nφ`` respectively.
Then the Cartesian componentes of each direction are obtained with

```math
[sin(θᵢ)*cos(φᵢ), sin(θᵢ)*sin(φᵢ), cos(θᵢ)]
```
Note that by discretization, the polar point may be computed more than once, in
general. Removing duplicates takes time; you can choose to remove duplicates or
not by setting the keyword argument `remove_duplicates` in the constructor
(default `true`).

### Examples

A `SphericalDirections` can be built in different ways. If you pass only one integer,
it is used to discretize both ``θ`` and ``φ``:

```jldoctest spherical_directions
julia> using LazySets.Approximations: SphericalDirections

julia> sd = SphericalDirections(3)
SphericalDirections(3, 3, Array{Float64,1}[[-1.22465e-16, -0.0, -1.0], [0.0, 0.0, 1.0]])

julia> sd.Nθ, sd.Nφ 
(3, 3)
```

Repeated directions are removed by default. Use the keyword argument `remove_duplicates`
to control this behavior:

```jldoctest spherical_directions
julia> sd_without_duplicates = SphericalDirections(3, remove_duplicates=true);

julia> length(sd_without_duplicates)
2

julia> sd_with_duplicates = SphericalDirections(3, remove_duplicates=false);

julia> length(sd_with_duplicates)
9
```

Pass two integers to control the discretization in ``θ`` and in ``φ`` separately:
```jldoctest spherical_directions
julia> sd_4_5 = SphericalDirections(4, 5);

julia> length(sd_4_5)
5

julia> sd_4_8 = SphericalDirections(4, 8);

julia> length(sd_4_8)
15
```
"""
struct SphericalDirections{N} <: AbstractDirections{N}
    Nθ::Int
    Nφ::Int
    stack::Vector{Vector{N}} # stores the spherical directions

    function SphericalDirections(Nθ::Int, Nφ::Int; N::Type{<:AbstractFloat}=Float64, remove_duplicates=true)
        if Nθ <= 1 || Nφ <= 1
            throw(ArgumentError("(Nθ, Nφ) = ($Nθ, $Nφ) is invalid; both shoud be at least 2"))
        end
        stack = Vector{Vector{N}}()
        θ = Compat.range(N(-pi), N(pi), length=Nθ)      # discretization of the azimuthal angle
        φ = Compat.range(N(0.0), N(2*pi), length=Nφ)    # discretization of the polar angle
        for φᵢ in φ
            for θⱼ in θ
                d = [sin(θⱼ)*cos(φᵢ), sin(θⱼ)*sin(φᵢ), cos(θⱼ)]
                if remove_duplicates && _isapprox_included(d, stack)
                    continue
                end
                push!(stack, d)
            end
        end
        return new{N}(Nθ, Nφ, stack)
    end
end

# convenience constructors
SphericalDirections(Nθ::Int; kwargs...) = SphericalDirections(Nθ, Nθ; kwargs...)

# common functions
Base.eltype(::Type{SphericalDirections{N}}) where {N<:AbstractFloat} = Vector{N}
Base.length(sd::SphericalDirections) = length(sd.stack)
dim(::SphericalDirections) = 3

@static if VERSION < v"0.7-"
    @eval begin
        Base.start(sd::SphericalDirections) = sd[1]
        Base.next(sd::SphericalDirections, state::Int) = (sd.stack[state], state+1)
        Base.done(sd::SphericalDirections, state) = state == length(sd.stack)+1
    end
else
    @eval begin
        function Base.iterate(sd::SphericalDirections, state::Int=1)
            state == length(sd.stack)+1 && return nothing
            return (sd.stack[state], state + 1)
        end
    end
end
