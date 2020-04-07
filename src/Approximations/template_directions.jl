import LazySets: dim, inner
using LazySets: isapproxzero

"""
    AbstractDirections{N, VN}

Abstract type for template direction representations.

### Notes

This type is parameterzed by `N` and `VN`, where:

- `N` stands for the numeric type
- `VN` stands for the vector type with coefficients of type `N`

Each subtype is an iterator over a set of prescribed directions.

All subtypes should implement the standard iterator methods from `Base`, namely
`Base.length` (returns the number of directions in the template), and `Base.iterate`.
Moreover, the following methods should be implemented:

- `dim`    -- return the ambient dimension of the template
- `eltype` -- return the type of each vector in the template

Optionally, subtypes may implement:

- `isbounding`   -- (defaults to `false`) return `true` if an overapproximation with
                    a list of template directions results in a bounded set, given a
                    bounded input set, and `false` otherwise
- `isnormalized` -- (defaults to `false`) returns `true` if each direction in the
                    given template has norm one w.r.t. the usual vector 2-norm
"""
abstract type AbstractDirections{N, VN} end

"""
    dim(ad::AbstractDirections)

Returns the dimension of the generated directions.

### Input

- `ad` -- template directions

### Output

The dimension of the generated directions.
"""
function dim(ad::AbstractDirections) end

"""
    isbounding(ad::Type{<:AbstractDirections})

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
function isbounding(ad::Type{<:AbstractDirections})
    return false
end

isbounding(::AD) where {AD<:AbstractDirections} = isbounding(AD)

"""
    isnormalized(ad::Type{<:AbstractDirections})

Returns whether the given template directions is normalized with respect to the
2-norm.

### Input

- `ad` -- template directions

### Output

`true` if the 2-norm of each element in `ad` is one and `false` otherwise.
"""
function isnormalized(ad::Type{<:AbstractDirections})
    return false
end

isnormalized(::AD) where {AD<:AbstractDirections} = isnormalized(AD)

# ==================================================
# Box directions
# ==================================================

"""
    BoxDirections{N, VN} <: AbstractDirections{N, VN}

Box directions representation.

### Fields

- `n` -- dimension

### Notes

Box directions can be seen as the vectors where only one entry is ±1, and all
other entries are 0. In dimension ``n``, there are ``2n`` such directions.

The deafault vector representation used in this template is a
`LazySets.Arrays.SingleEntryVector`, although other implementations can be used
such as a regular `Vector` and a sparse vector, `SparseVector`.`

### Examples

The template can be constructed by passing the dimension. For example, in dimension
two,

```jldoctest dirs_Box
julia> dirs = BoxDirections(2)
BoxDirections{Float64,LazySets.Arrays.SingleEntryVector{Float64}}(2)

julia> length(dirs)
4
```

By default, each direction is represented in this iterator as a `SingleEntryVector`,
i.e. a vector with only one non-zero element,

```jldoctest dirs_Box
julia> eltype(dirs)
LazySets.Arrays.SingleEntryVector{Float64}
```
In two dimensions, the directions defined by `BoxDirections` are normal to the
facets of a box.

```jldoctest dirs_Box
julia> collect(dirs)
4-element Array{LazySets.Arrays.SingleEntryVector{Float64},1}:
 [1.0, 0.0]
 [0.0, 1.0]
 [0.0, -1.0]
 [-1.0, 0.0]
```

The numeric type can be specified as well:

```jldoctest
julia> BoxDirections{Rational{Int}}(10)
BoxDirections{Rational{Int64},LazySets.Arrays.SingleEntryVector{Rational{Int64}}}(10)

julia> length(ans)
20
```
"""
struct BoxDirections{N, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}
    n::Int
end

# convenience constructor for type Float64
BoxDirections(n::Int) = BoxDirections{Float64, SingleEntryVector{Float64}}(n)

# constructor where only N is specified
BoxDirections{N}(n::Int) where {N} = BoxDirections{N, SingleEntryVector{N}}(n)

Base.eltype(::Type{BoxDirections{N, VN}}) where {N, VN} = VN
Base.length(bd::BoxDirections) = 2 * bd.n

# interface functions
dim(bd::BoxDirections) = bd.n
isbounding(::Type{<:BoxDirections}) = true
isnormalized(::Type{<:BoxDirections}) = true

# The idea is that a positive states run through vectors wtih +1 entry,
# and negative states run through -1 entries
# (1, 0)   state = 1
# (0, 1)   state = 2
# (0, -1)  state = -2
# (-1, 0)  state = -1
function Base.iterate(bd::BoxDirections{N, SingleEntryVector{N}}, state::Int=1) where {N}
    if state == 0
        return nothing
    end
    vec = SingleEntryVector(abs(state), bd.n, convert(N, sign(state)))
    state = (state == bd.n) ? -bd.n : state + 1
    return (vec, state)
end

function Base.iterate(bd::BoxDirections{N, Vector{N}}, state::Int=1) where {N}
    if state == 0
        return nothing
    end
    vec = zeros(N, bd.n)
    vec[abs(state)] = convert(N, sign(state))
    state = (state == bd.n) ? -bd.n : state + 1
    return (vec, state)
end

function Base.iterate(bd::BoxDirections{N, SparseVector{N, Int}}, state::Int=1) where {N}
    if state == 0
        return nothing
    end
    vec = sparsevec([abs(state)], convert(N, sign(state)), bd.n)
    state = (state == bd.n) ? -bd.n : state + 1
    return (vec, state)
end

# ==================================================
# Octagonal directions
# ==================================================

"""
    OctDirections{N, VN} <: AbstractDirections{N, VN}

Octagon directions representation.

### Fields

- `n` -- dimension

### Notes

Octagon directions consist of all vectors that are zero almost everywhere except
in two dimensions ``i``, ``j`` (possibly ``i = j``) where it is ``±1``. In dimension
``n``, there are ``2n^2`` such directions.

## Examples

The template can be constructed by passing the dimension. For example, in dimension
two,

```jldoctest dirs_Oct
julia> dirs = OctDirections(2)
OctDirections{Float64,SparseVector{Float64,Int64}}(2)

julia> length(dirs) # number of directions
8
```
By default, each direction is represented in this iterator as a sparse vector:

```jldoctest dirs_Oct
julia> eltype(dirs)
SparseVector{Float64,Int64}
```
In two dimensions, the directions defined by `BoxDiagDirections` are normal to
the facets of an octagon.

```jldoctest dirs_Oct
julia> first(dirs)
2-element SparseVector{Float64,Int64} with 2 stored entries:
  [1]  =  1.0
  [2]  =  1.0

julia> Vector.(collect(dirs))
8-element Array{Array{Float64,1},1}:
 [1.0, 1.0]
 [1.0, -1.0]
 [-1.0, 1.0]
 [-1.0, -1.0]
 [1.0, 0.0]
 [0.0, 1.0]
 [0.0, -1.0]
 [-1.0, 0.0]
```

The numeric type can be specified as well:

```jldoctest
julia> OctDirections{Rational{Int}}(10)
OctDirections{Rational{Int64},SparseVector{Rational{Int64},Int64}}(10)

julia> length(ans)
200
```
"""
struct OctDirections{N, VN} <: AbstractDirections{N, VN}
    n::Int
end

# constructor for type Float64
OctDirections(n::Int) = OctDirections{Float64, SparseVector{Float64, Int}}(n)

# constructor where only N is specified
OctDirections{N}(n::Int) where {N} = OctDirections{N, SparseVector{N, Int}}(n)

Base.eltype(::Type{OctDirections{N, VN}}) where {N, VN} = VN
Base.length(od::OctDirections) = 2 * od.n^2

# interface functions
dim(od::OctDirections) = od.n
isbounding(::Type{<:OctDirections}) = true
isnormalized(::Type{<:OctDirections}) = false

function Base.iterate(od::OctDirections{N, SparseVector{N, Int}}) where {N}
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
function Base.iterate(od::OctDirections{N, SparseVector{N, Int}}, state::Tuple) where {N}
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

function Base.iterate(od::OctDirections{N, SparseVector{N, Int}}, state::Int) where {N}
    # continue with box directions
    return iterate(BoxDirections{N, SparseVector{N, Int}}(od.n), state)
end

# ==================================================
# Box-diagonal directions
# ==================================================

"""
    BoxDiagDirections{N, VN} <: AbstractDirections{N, VN}

Box-diagonal directions representation.

### Fields

- `n` -- dimension

### Notes

Box-diagonal directions can be seen as the union of diagonal directions (all
entries are ±1) and box directions (one entry is ±1, all other entries are 0).
The iterator first enumerates all diagonal directions, and then all box
directions. In dimension ``n``, there are in total ``2^n + 2n`` such directions.

## Examples

The template can be constructed by passing the dimension. For example, in dimension
two,

```jldoctest dirs_BoxDiag
julia> dirs = BoxDiagDirections(2)
BoxDiagDirections{Float64,Array{Float64,1}}(2)

julia> length(dirs) # number of directions
8
```
By default, each direction is represented in this iterator as a regular vector:

```jldoctest dirs_BoxDiag
julia> eltype(dirs)
Array{Float64,1}
```
In two dimensions, the directions defined by `BoxDiagDirections` are normal to
the facets of an octagon.

```jldoctest dirs_BoxDiag
julia> collect(dirs)
8-element Array{Array{Float64,1},1}:
 [1.0, 1.0]
 [-1.0, 1.0]
 [1.0, -1.0]
 [-1.0, -1.0]
 [1.0, 0.0]
 [0.0, 1.0]
 [0.0, -1.0]
 [-1.0, 0.0]
```

The numeric type can be specified as well:

```jldoctest
julia> BoxDiagDirections{Rational{Int}}(10)
BoxDiagDirections{Rational{Int64},Array{Rational{Int64},1}}(10)

julia> length(ans)
1044
```
"""
struct BoxDiagDirections{N, VN} <: AbstractDirections{N, VN}
    n::Int
end

# constructor for type Float64
BoxDiagDirections(n::Int) = BoxDiagDirections{Float64, Vector{Float64}}(n)

# constructor where only N is specified
BoxDiagDirections{N}(n::Int) where {N} = BoxDiagDirections{N, Vector{N}}(n)

Base.eltype(::Type{BoxDiagDirections{N, VN}}) where {N, VN} = VN
Base.length(bdd::BoxDiagDirections) = bdd.n == 1 ? 2 : 2^bdd.n + 2 * bdd.n

# interface function
dim(bdd::BoxDiagDirections) = bdd.n
isbounding(::Type{<:BoxDiagDirections}) = true
isnormalized(::Type{<:BoxDiagDirections}) = false

function Base.iterate(bdd::BoxDiagDirections{N, Vector{N}}) where {N}
    return (ones(N, bdd.n), ones(N, bdd.n))
end

function Base.iterate(bdd::BoxDiagDirections{N}, state::Vector{N}) where {N}
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

function Base.iterate(bdd::BoxDiagDirections{N, Vector{N}}, state::Int) where {N}
    # continue with box directions
    return iterate(BoxDirections{N, Vector{N}}(bdd.n), state)
end

# ==================================================
# Polar directions
# ==================================================

"""
    PolarDirections{N<:AbstractFloat, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}

Polar directions representation.

### Fields

- `Nφ` -- length of the partition of the polar angle

### Notes

The `PolarDirections` constructor provides a sample of the unit sphere
in ``\\mathbb{R}^2``, which is parameterized by the polar angle
``φ ∈ Dφ := [0, 2π]``; see the wikipedia entry
[Polar coordinate system](https://en.wikipedia.org/wiki/Polar_coordinate_system)
for details.

The integer argument ``Nφ`` defines how many samples of ``Dφ`` are taken. The
Cartesian components of each direction are obtained with

```math
[cos(φᵢ), sin(φᵢ)].
```

### Examples

The integer passed as an argument is used to discretize ``φ``:

```jldoctest; filter = r"2246[0-9]*e-16"
julia> pd = PolarDirections(2)
PolarDirections{Float64,Array{Float64,1}}(2, [[1.0, 0.0], [-1.0, 1.2246467991473532e-16]])

julia> length(pd)
2
```
"""
struct PolarDirections{N<:AbstractFloat, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}
    Nφ::Int
    stack::Vector{VN} # stores the polar directions
end

# convenience constructor for type Float64
PolarDirections(Nφ::Int) = PolarDirections{Float64, Vector{Float64}}(Nφ)

function PolarDirections{Float64, Vector{Float64}}(Nφ::Int)
    if Nφ <= 0
        throw(ArgumentError("Nφ = $Nφ is invalid; it shoud be at least 1"))
    end
    stack = Vector{Vector{Float64}}(undef, Nφ)
    φ = range(0.0, 2*pi, length=Nφ+1)  # discretization of the polar angle

    @inbounds for i in 1:Nφ  # skip last (repeated) angle
        stack[i] = [cos(φ[i]), sin(φ[i])]
    end
    return PolarDirections{Float64, Vector{Float64}}(Nφ, stack)
end

Base.eltype(::Type{PolarDirections{N, VN}}) where {N, VN} = VN
Base.length(pd::PolarDirections) = pd.Nφ

# interface functions
dim(pd::PolarDirections) = 2
isbounding(pd::Type{<:PolarDirections}) = pd.Nφ > 2
isnormalized(::Type{<:PolarDirections}) = true

function Base.iterate(pd::PolarDirections{N, Vector{N}}, state::Int=1) where {N}
    state == pd.Nφ + 1 && return nothing
    return (pd.stack[state], state + 1)
end

# ==================================================
# Spherical directions
# ==================================================

"""
    SphericalDirections{N<:AbstractFloat, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}

Spherical directions representation.

### Fields

- `Nθ`    -- length of the partition of the azimuthal angle
- `Nφ`    -- length of the partition of the polar angle
- `stack` -- list of computed directions

### Notes

The `SphericalDirections` constructor provides a sample of the unit sphere
in ``\\mathbb{R}^3``, which is parameterized by the azimuthal and polar angles
``θ ∈ Dθ := [0, π]`` and ``φ ∈ Dφ := [0, 2π]`` respectively, see the wikipedia
entry [Spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
for details.

The integer arguments ``Nθ`` and ``Nφ`` define how many samples along the domains
``Dθ`` and ``Dφ`` respectively are taken. The Cartesian components of each direction
are obtained with

```math
[sin(θᵢ)*cos(φᵢ), sin(θᵢ)*sin(φᵢ), cos(θᵢ)].
```

The north and south poles are treated separately so that those points
are not considered more than once.

### Examples

A `SphericalDirections` template can be built in different ways. If you pass
only one integer, the same value is used to discretize both ``θ`` and ``φ``:

```jldoctest spherical_directions; filter = r"1232[0-9]*e-17.*2246[0-9]*e-16.*1232[0-9]*e-17"
julia> sd = SphericalDirections(3)
SphericalDirections{Float64,Array{Float64,1}}(3, 3, [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [1.0, 0.0, 6.123233995736766e-17], [-1.0, 1.2246467991473532e-16, 6.123233995736766e-17]])

julia> sd.Nθ, sd.Nφ
(3, 3)

julia> length(sd)
4
```

Pass two integers to control the discretization in ``θ`` and in ``φ`` separately:

```jldoctest spherical_directions
julia> sd = SphericalDirections(4, 5);

julia> length(sd)
10

julia> sd = SphericalDirections(4, 8);

julia> length(sd)
16
```
"""
struct SphericalDirections{N<:AbstractFloat, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}
    Nθ::Int
    Nφ::Int
    stack::Vector{VN} # stores the spherical directions
end

# convenience constructors
SphericalDirections(Nθ::Int) = SphericalDirections(Nθ, Nθ)
SphericalDirections(Nθ::Int, Nφ::Int) = SphericalDirections{Float64, Vector{Float64}}(Nθ::Int, Nφ::Int)

function SphericalDirections{Float64, Vector{Float64}}(Nθ::Int, Nφ::Int)
    if Nθ <= 1 || Nφ <= 1
        throw(ArgumentError("(Nθ, Nφ) = ($Nθ, $Nφ) is invalid; both shoud be at least 2"))
    end
    stack = Vector{Vector{Float64}}()
    θ = range(0.0, pi, length=Nθ)    # discretization of the azimuthal angle
    φ = range(0.0, 2*pi, length=Nφ)  # discretization of the polar angle

    # add north pole (θ = 0)
    push!(stack, Float64[0, 0, 1])

    # add south pole (θ = pi)
    push!(stack, Float64[0, 0, -1])

    for φᵢ in φ[1:Nφ-1]  # delete repeated angle
        for θⱼ in θ[2:Nθ-1] # delete north and south poles
            d = [sin(θⱼ)*cos(φᵢ), sin(θⱼ)*sin(φᵢ), cos(θⱼ)]
            push!(stack, d)
        end
    end
    return SphericalDirections{Float64, Vector{Float64}}(Nθ, Nφ, stack)
end

Base.eltype(::Type{SphericalDirections{N, VN}}) where {N, VN} = VN
Base.length(sd::SphericalDirections) = length(sd.stack)

# interface functions
dim(::SphericalDirections) = 3
isbounding(sd::Type{<:SphericalDirections}) = sd.Nθ > 2 && sd.Nφ > 2
isnormalized(::Type{<:SphericalDirections}) = true

function Base.iterate(sd::SphericalDirections, state::Int=1)
    state == length(sd.stack) + 1 && return nothing
    return (sd.stack[state], state + 1)
end

# ==================================================
# Custom template directions
# ==================================================

"""
    CustomDirections{N, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}

User-defined template directions.

### Fields

- `directions` -- list of template directions
- `n`          -- dimension
- `isbounding` -- boundedness status

### Notes

This struct is a wrapper type for a set of user-defined directions which are
iterated over. It has fields for the list of directions, the set dimension,
and (boolean) cache fields for the boundedness and normalization properties.
The latter are checked by default upon construction.

To check boundedness, we overapproximate the unit ball in the infinity
norm using the given directions and check if the resulting set is bounded.

The dimension will also be determined automatically, unless the empty vector is
passed (in which case the optional argument `n` needs to be specified).

## Examples

Creating a template with box directions in dimension two:

```jldoctest
julia> dirs = CustomDirections([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]])
CustomDirections{Float64,Array{Float64,1}}([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]], 2, true, true)

julia> LazySets.Approximations.isbounding(dirs)
true

julia> LazySets.Approximations.isnormalized(dirs)
true
```
"""
struct CustomDirections{N, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}
    directions::Vector{VN}
    n::Int
    bounded::Bool
    normalized::Bool
end

function CustomDirections(directions::Vector{VN},
                          n::Int=_determine_dimension(directions);
                          check_boundedness::Bool=true,
                          check_normalization::Bool=true) where {N, VN<:AbstractVector{N}}
    bounded = check_boundedness ? _isbounding(directions) : false
    normalized = check_normalization ? _isnormalized(directions) : false
    return CustomDirections(directions, n, bounded, normalized)
end

function _determine_dimension(directions)
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

function _isnormalized(directions::Vector{VN}) where {N, VN<:AbstractVector{N}}
    return all(x -> _isapprox(x, one(N)), norm.(directions, 2))
end

Base.eltype(::Type{CustomDirections{N, VN}}) where {N, VN} = VN
Base.length(cd::CustomDirections) = length(cd.directions)

# interface functions
dim(cd::CustomDirections) = cd.n
isbounding(cd::Type{<:CustomDirections}) = false # it is not a property of the type
isbounding(cd::CustomDirections) = cd.bounded
isnormalized(cd::Type{<:CustomDirections}) = false # it is not a property of the type
isnormalized(cd::CustomDirections) = cd.normalized

function Base.iterate(cd::CustomDirections{N}, state::Int=1) where {N}
    if state > length(cd.directions)
        return nothing
    end
    vec = cd.directions[state]
    state = state + 1
    return (vec, state)
end
