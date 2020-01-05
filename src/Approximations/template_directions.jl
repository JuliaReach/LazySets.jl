import LazySets: dim, inner
using LazySets: isapproxzero

"""
    AbstractDirections{N}

Abstract type for template direction representations.

### Notes

All subtypes should implement the standard iterator methods from `Base` and the
function `dim(d<:AbstractDirections)`.
"""
abstract type AbstractDirections{N} end

"""
    dim(ad::AbstractDirections)

Returns the dimension of the generated directions.

### Input

- `ad` -- template directions

### Output

The dimension of the generated directions.
"""
function dim(ad::AbstractDirections)
    return ad.n
end

"""
    isbounding(ad::AbstractDirections)

Checks if an overapproximation with a list of template directions results in a
bounded set, given a bounded input set.

### Input

- `ad` -- template directions

### Output

Given a bounded set ``X``, we can construct an outer approximation of ``X`` by
using the template directions `ad` as normal vectors of the facets.
If this function returns `true`, then the result is again a bounded set (i.e., a
polytope).
Note that the result does not depend on the specific shape of ``X``, as long as
``X`` is bounded.

### Notes

By default, this function returns `false` in order to be conservative.
Custom subtypes of `AbstractDirections` should hence add a method for this
function.
"""
function isbounding(ad::AbstractDirections)
    return false
end

# ==================================================
# Box directions
# ==================================================

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
isbounding(::BoxDirections) = true

function Base.iterate(bd::BoxDirections{N}, state::Int=1) where {N}
    if state == 0
        return nothing
    end
    vec = SingleEntryVector(abs(state), bd.n, convert(N, sign(state)))
    state = (state == bd.n) ? -bd.n : state + 1
    return (vec, state)
end

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
isbounding(::OctDirections) = true

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
isbounding(::BoxDiagDirections) = true

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

# ==================================================
# Polar directions
# ==================================================

"""
    PolarDirections{N<:AbstractFloat} <: AbstractDirections{N}

Polar directions representation.

### Fields

- `Nφ`    -- length of the partition of the polar angle

### Notes

The `PolarDirections` constructor provides a sample of the unit sphere
in ``\\mathbb{R}^2``, which is parameterized by the polar angles
``φ ∈ Dφ := [0, 2π]`` respectively; see the wikipedia entry
[Polar coordinate system](https://en.wikipedia.org/wiki/Polar_coordinate_system).
The domain ``Dφ`` is discretized in ``Nφ`` pieces.
Then the Cartesian components of each direction are obtained with

```math
[cos(φᵢ), sin(φᵢ)].
```

### Examples

The integer passed as an argument is used to discretize ``φ``:

```jldoctest; filter = r"2246[0-9]*e-16"
julia> using LazySets.Approximations: PolarDirections

julia> pd = PolarDirections(2)
PolarDirections{Float64}(2, Array{Float64,1}[[1.0, 0.0], [-1.0, 1.22465e-16]])

julia> pd.Nφ
2
```
"""
struct PolarDirections{N<:AbstractFloat} <: AbstractDirections{N}
    Nφ::Int
    stack::Vector{Vector{N}} # stores the polar directions

    function PolarDirections{N}(Nφ::Int) where {N<:AbstractFloat}
        if Nφ <= 0
            throw(ArgumentError("Nφ = $Nφ is invalid; it shoud be at least 1"))
        end
        stack = Vector{Vector{N}}()
        # discretization of the polar angle
        φ = range(N(0.0), N(2*pi), length=Nφ+1)

        for φᵢ in φ[1:Nφ]  # skip last (repeated) angle
            d = N[cos(φᵢ), sin(φᵢ)]
            push!(stack, d)
        end
        return new{N}(Nφ, stack)
    end
end

PolarDirections(Nφ::Int) = PolarDirections{Float64}(Nφ)

# common functions
Base.eltype(::Type{PolarDirections{N}}) where {N} = Vector{N}
Base.length(pd::PolarDirections) = length(pd.stack)
dim(::PolarDirections) = 2
isbounding(pd::PolarDirections) = pd.Nφ > 2

function Base.iterate(pd::PolarDirections, state::Int=1)
    state == length(pd.stack)+1 && return nothing
    return (pd.stack[state], state + 1)
end

# ==================================================
# Spherical directions
# ==================================================

"""
    SphericalDirections{N<:AbstractFloat} <: AbstractDirections{N}

Spherical directions representation.

### Fields

- `Nθ`    -- length of the partition of the azimuthal angle
- `Nφ`    -- length of the partition of the polar angle
- `stack` -- list of computed directions

### Notes

The `SphericalDirections` constructor provides a sample of the unit sphere
in ``\\mathbb{R}^3``, which is parameterized by the azimuthal and polar angles
``θ ∈ Dθ := [0, π]`` and ``φ ∈ Dφ := [0, 2π]`` respectively, see the wikipedia
entry [Spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system).
The domains ``Dθ`` and ``Dφ`` are discretized in ``Nθ`` and ``Nφ`` respectively.
Then the Cartesian components of each direction are obtained with

```math
[sin(θᵢ)*cos(φᵢ), sin(θᵢ)*sin(φᵢ), cos(θᵢ)].
```

The north and south poles are treated separately so that those points
are not considered more than once.

### Examples

A `SphericalDirections` can be built in different ways. If you pass only one integer,
it is used to discretize both ``θ`` and ``φ``:

```jldoctest spherical_directions; filter = r"1232[0-9]*e-17.*2246[0-9]*e-16.*1232[0-9]*e-17"
julia> using LazySets.Approximations: SphericalDirections

julia> sd = SphericalDirections(3)
SphericalDirections{Float64}(3, 3, Array{Float64,1}[[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [1.0, 0.0, 6.12323e-17], [-1.0, 1.22465e-16, 6.12323e-17]])

julia> sd.Nθ, sd.Nφ 
(3, 3)
```

Pass two integers to control the discretization in ``θ`` and in ``φ`` separately:

```jldoctest spherical_directions
julia> sd_4_5 = SphericalDirections(4, 5);

julia> length(sd_4_5)
10

julia> sd_4_8 = SphericalDirections(4, 8);

julia> length(sd_4_8)
16
```
"""
struct SphericalDirections{N<:AbstractFloat} <: AbstractDirections{N}
    Nθ::Int
    Nφ::Int
    stack::Vector{Vector{N}} # stores the spherical directions

    function SphericalDirections{N}(Nθ::Int, Nφ::Int) where {N<:AbstractFloat}
        if Nθ <= 1 || Nφ <= 1
            throw(ArgumentError("(Nθ, Nφ) = ($Nθ, $Nφ) is invalid; both shoud be at least 2"))
        end
        stack = Vector{Vector{N}}()
        θ = range(N(0), N(pi), length=Nθ)      # discretization of the azimuthal angle
        φ = range(N(0.0), N(2*pi), length=Nφ)  # discretization of the polar angle

        # add north pole (θ = 0)
        push!(stack, N[0, 0, 1])

        # add south pole (θ = pi)
        push!(stack, N[0, 0, -1])

        for φᵢ in φ[1:Nφ-1]  # delete repeated angle
            for θⱼ in θ[2:Nθ-1] # delete north and south poles
                d = N[sin(θⱼ)*cos(φᵢ), sin(θⱼ)*sin(φᵢ), cos(θⱼ)]
                push!(stack, d)
            end
        end
        return new{N}(Nθ, Nφ, stack)
    end
end

# convenience constructors
SphericalDirections(Nθ::Int) = SphericalDirections{Float64}(Nθ, Nθ)
SphericalDirections(Nθ::Int, Nφ::Int) = SphericalDirections{Float64}(Nθ, Nφ)

# common functions
Base.eltype(::Type{SphericalDirections{N}}) where {N} = Vector{N}
Base.length(sd::SphericalDirections) = length(sd.stack)
dim(::SphericalDirections) = 3
isbounding(sd::SphericalDirections) = sd.Nθ > 2 && sd.Nφ > 2

function Base.iterate(sd::SphericalDirections, state::Int=1)
    state == length(sd.stack)+1 && return nothing
    return (sd.stack[state], state + 1)
end

# ==================================================
# Custom template directions
# ==================================================

"""
    CustomDirections{N, VN<:AbstractVector{N}} <: AbstractDirections{N}

User-defined template directions.

### Fields

- `directions` -- list of template directions
- `n`          -- dimension
- `isbounding` -- boundedness status

### Notes

Since custom directions may not be bounded, boundedness can either be asserted
in the constructor or it will be determined automatically by performing an
experiment (more precisely, we overapproximate the unit ball in the infinity
norm using the given directions).

The dimension will also be determined automatically, unless the empty vector is
passed (in which case the optional argument `n` needs to be specified).
"""
struct CustomDirections{N, VN<:AbstractVector{N}} <: AbstractDirections{N}
    directions::Vector{VN}
    n::Int
    isbounding::Bool

    CustomDirections{N, VN}(directions::Vector{VN};
                            n::Int=determine_dimension(directions),
                            isbounding::Bool=_isbounding(directions)
                           ) where {N, VN<:AbstractVector{N}} =
        new{N, VN}(directions, n, isbounding)
end

# convenience constructor
CustomDirections(directions::Vector{VN};
                 n::Int=determine_dimension(directions),
                 isbounding::Bool=_isbounding(directions)
                ) where {N, VN<:AbstractVector{N}} =
    CustomDirections{N, VN}(directions; n=n, isbounding=isbounding)

function determine_dimension(directions)
    isempty(directions) && throw(ArgumentError("empty template directions " *
                                               "need a specified dimension"))
    return length(directions[1])
end

function _isbounding(directions::Vector{VN}) where {N, VN<:AbstractVector{N}}
    isempty(directions) && return false

    # check boundedness of the polyhedron `⋂_a ax <= 1` where `a` is a direction
    P = HPolyhedron([HalfSpace(dir, one(N)) for dir in directions])
    return isbounded(P)
end

Base.eltype(::Type{CustomDirections{N, VN}}) where {N, VN} = VN
Base.length(cd::CustomDirections) = length(cd.directions)
isbounding(cd::CustomDirections) = cd.isbounding

function Base.iterate(cd::CustomDirections{N}, state::Int=1) where {N}
    if state > length(cd.directions)
        return nothing
    end
    vec = cd.directions[state]
    state = state + 1
    return (vec, state)
end
