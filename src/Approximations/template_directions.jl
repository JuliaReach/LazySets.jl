import LazySets: dim, inner
using LazySets: isapproxzero

"""
    AbstractDirections{N, VN}

Abstract type for representations of direction vectors.

### Notes

This type is parameterized by `N` and `VN`, where:

- `N` stands for the numeric type
- `VN` stands for the vector type with coefficients of type `N`

Each implementing subtype is an iterator over a set of directions.
For that they implement the standard iterator methods from `Base`, namely
`Base.length` (returns the number of directions) and `Base.iterate`.
Moreover, the following methods should be implemented:

- `dim`    -- return the ambient dimension of the vectors
- `eltype` -- return the type of each vector

Optionally, subtypes may implement:

- `isbounding`   -- (defaults to `false`) return `true` if an overapproximation
                    with the direction vectors results in a bounded set, given a
                    bounded input set, and `false` otherwise
- `isnormalized` -- (defaults to `false`) is `true` if each direction vector has
                    norm one w.r.t. the usual vector 2-norm
"""
abstract type AbstractDirections{N, VN} end

"""
    dim(ad::AbstractDirections)

Returns the dimension of the generated directions.

### Input

- `ad` -- direction vectors

### Output

The ambient dimension of the generated directions.
"""
function dim(ad::AbstractDirections) end

"""
    isbounding(ad::AbstractDirections)
    isbounding(ad::Type{<:AbstractDirections})

Check whether an overapproximation with a set of direction vectors results in a
bounded set, given a bounded input set.

### Input

- `ad` -- direction vectors or a subtype of `AbstractDirections`

### Output

Given a bounded set ``X``, we can construct an outer polyhedral approximation of
``X`` by using the direction vectors `ad` as normal vectors of the facets.
If this function returns `true`, then the result is again guaranteed to be a
bounded set (i.e., a polytope).
Note that the result does not depend on the specific shape of ``X``, as long as
``X`` is bounded.

### Notes

By default, this function returns `false` in order to be conservative.
Custom subtypes of `AbstractDirections` should hence add a method for this
function.

The function can be applied to an instance of an `AbstractDirections` subtype or
to the subtype itself. By default, the check on the instance falls back to the
check on the subtype.
"""
function isbounding(ad::Type{<:AbstractDirections})
    return false
end

isbounding(::AD) where {AD<:AbstractDirections} = isbounding(AD)

"""
    isnormalized(ad::AbstractDirections)
    isnormalized(ad::Type{<:AbstractDirections})

Check whether the given direction vectors are normalized with respect to the
2-norm.

### Input

- `ad` -- direction vectors or a subtype of `AbstractDirections`

### Output

`true` if the 2-norm of each element in `ad` is one and `false` otherwise.

### Notes

By default, this function returns `false` in order to be conservative.
Custom subtypes of `AbstractDirections` should hence add a method for this
function.

The function can be applied to an instance of an `AbstractDirections` subtype or
to the subtype itself. By default, the check on the instance falls back to the
check on the subtype.
"""
function isnormalized(ad::Type{<:AbstractDirections})
    return false
end

isnormalized(::AD) where {AD<:AbstractDirections} = isnormalized(AD)

"""
    project(S::LazySet,
            block::AbstractVector{Int},
            directions::Type{<:AbstractDirections},
            [n]::Int;
            [kwargs...]
           )

Project a high-dimensional set to a given block using direction vectors.

### Input

- `S`          -- set
- `block`      -- block structure - a vector with the dimensions of interest
- `directions` -- direction vectors
- `n`          -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

The polyhedral overapproximation of the projection of `S` in the given
directions.
"""
@inline function project(S::LazySet,
                         block::AbstractVector{Int},
                         directions::Type{<:AbstractDirections},
                         n::Int=dim(S);
                         kwargs...
                        )
    lm = project(S, block, LinearMap, n; kwargs...)
    return overapproximate(lm, directions(length(block)))
end

"""
    BoxDirections{N, VN} <: AbstractDirections{N, VN}

Box directions representation.

### Fields

- `n` -- dimension

### Notes

Box directions can be seen as the vectors where only one entry is ±1, and all
other entries are 0. In dimension ``n``, there are ``2n`` such directions.

The default vector representation used in this template is a
`ReachabilityBase.Arrays.SingleEntryVector`, although other implementations can
be used such as a regular `Vector` and a `SparseVector`.

### Examples

The template can be constructed by passing the dimension. For example, in
dimension two:

```jldoctest dirs_Box
julia> dirs = BoxDirections(2)
BoxDirections{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}(2)

julia> length(dirs)
4
```

By default, each direction is represented as a `SingleEntryVector`, i.e., a
vector with only one non-zero element,

```jldoctest dirs_Box
julia> eltype(dirs)
ReachabilityBase.Arrays.SingleEntryVector{Float64}
```
In two dimensions, the directions defined by `BoxDirections` are normal to the
facets of a box.

```jldoctest dirs_Box
julia> collect(dirs)
4-element Vector{ReachabilityBase.Arrays.SingleEntryVector{Float64}}:
 [1.0, 0.0]
 [0.0, 1.0]
 [0.0, -1.0]
 [-1.0, 0.0]
```

The numeric type can be specified as well:

```jldoctest
julia> BoxDirections{Rational{Int}}(10)
BoxDirections{Rational{Int64}, ReachabilityBase.Arrays.SingleEntryVector{Rational{Int64}}}(10)

julia> length(ans)
20
```
"""
struct BoxDirections{N, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}
    n::Int
end

# constructor where only N is specified
BoxDirections{N}(n::Int) where {N} = BoxDirections{N, SingleEntryVector{N}}(n)

# convenience constructor for type Float64
BoxDirections(n::Int) = BoxDirections{Float64}(n)

Base.eltype(::Type{BoxDirections{N, VN}}) where {N, VN} = VN
Base.length(bd::BoxDirections) = 2 * bd.n

# interface functions
dim(bd::BoxDirections) = bd.n
isbounding(::Type{<:BoxDirections}) = true
isnormalized(::Type{<:BoxDirections}) = true

# The idea is that positive states run through vectors with a +1 entry,
# and negative states run through vectors with a -1 entry
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
    @inbounds vec[abs(state)] = convert(N, sign(state))
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

"""
    OctDirections{N, VN} <: AbstractDirections{N, VN}

Octagon directions representation.

### Fields

- `n` -- dimension

### Notes

Octagon directions consist of all vectors that are zero almost everywhere except
in two dimensions ``i``, ``j`` (possibly ``i = j``) where it is ``±1``. In
dimension ``n``, there are ``2n^2`` such directions.

### Examples

The template can be constructed by passing the dimension. For example, in
dimension two:

```jldoctest dirs_Oct
julia> dirs = OctDirections(2)
OctDirections{Float64, SparseArrays.SparseVector{Float64, Int64}}(2)

julia> length(dirs) # number of directions
8
```
By default, the directions are represented as sparse vectors:

```jldoctest dirs_Oct
julia> eltype(dirs)
SparseArrays.SparseVector{Float64, Int64}
```
In two dimensions, the directions are normal to the facets of an octagon.

```jldoctest dirs_Oct
julia> first(dirs)
2-element SparseArrays.SparseVector{Float64, Int64} with 2 stored entries:
  [1]  =  1.0
  [2]  =  1.0

julia> Vector.(collect(dirs))
8-element Vector{Vector{Float64}}:
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
OctDirections{Rational{Int64}, SparseArrays.SparseVector{Rational{Int64}, Int64}}(10)

julia> length(ans)
200
```
"""
struct OctDirections{N, VN} <: AbstractDirections{N, VN}
    n::Int
end

# constructor where only N is specified
OctDirections{N}(n::Int) where {N} = OctDirections{N, SparseVector{N, Int}}(n)

# constructor for type Float64
OctDirections(n::Int) = OctDirections{Float64}(n)

Base.eltype(::Type{OctDirections{N, VN}}) where {N, VN} = VN
Base.length(od::OctDirections) = 2 * od.n^2

# interface functions
dim(od::OctDirections) = od.n
isbounding(::Type{<:OctDirections}) = true
isnormalized(::Type{<:OctDirections}) = false

function _zeros_oct(n, ::Type{<:SparseVector{N}}) where {N}
    return spzeros(N, n)
end

function Base.iterate(od::OctDirections{N, VN}) where {N, VN}
    if od.n == 1
        # fall back to box directions in 1D case
        return iterate(od, 1)
    end
    vec = _zeros_oct(od.n, VN)
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
function _iterate_state(od::OctDirections{N}, state) where {N}
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

function Base.iterate(od::OctDirections, state::Tuple)
    _iterate_state(od, state)
end

function Base.iterate(od::OctDirections{N, VN}, state::Int) where {N, VN}
    # continue with box directions
    return iterate(BoxDirections{N, VN}(od.n), state)
end

# implementation with regular Vector
function _zeros_oct(n, ::Type{Vector{N}}) where {N}
    return zeros(N, n)
end

"""
    DiagDirections{N, VN} <: AbstractDirections{N, VN}

Diagonal directions representation.

### Fields

- `n` -- dimension

### Notes

Diagonal directions are vectors where all entries are ±1. In dimension ``n``,
there are in total ``2^n`` such directions.

## Examples

The template can be constructed by passing the dimension. For example, in
dimension two:

```jldoctest dirs_Diag
julia> dirs = DiagDirections(2)
DiagDirections{Float64, Vector{Float64}}(2)

julia> length(dirs) # number of directions
4
```
By default, each direction is represented as a regular `Vector`:

```jldoctest dirs_Diag
julia> eltype(dirs)
Vector{Float64} (alias for Array{Float64, 1})
```
In two dimensions, the directions defined by `DiagDirections` are normal to the
facets of a ball in the 1-norm.

```jldoctest dirs_Diag
julia> collect(dirs)
4-element Vector{Vector{Float64}}:
 [1.0, 1.0]
 [-1.0, 1.0]
 [1.0, -1.0]
 [-1.0, -1.0]
```

The numeric type can be specified as well:

```jldoctest
julia> DiagDirections{Rational{Int}}(10)
DiagDirections{Rational{Int64}, Vector{Rational{Int64}}}(10)

julia> length(ans)
1024
```
"""
struct DiagDirections{N, VN} <: AbstractDirections{N, VN}
    n::Int
end

# constructor where only N is specified
DiagDirections{N}(n::Int) where {N} = DiagDirections{N, Vector{N}}(n)

# constructor for type Float64
DiagDirections(n::Int) = DiagDirections{Float64}(n)

Base.eltype(::Type{DiagDirections{N, VN}}) where {N, VN} = VN
Base.length(dd::DiagDirections) = 2^dd.n

# interface function
dim(dd::DiagDirections) = dd.n
isbounding(::Type{<:DiagDirections}) = true
isnormalized(::Type{<:DiagDirections}) = false

function Base.iterate(dd::DiagDirections{N, Vector{N}}) where {N}
    return (ones(N, dd.n), ones(N, dd.n))
end

function Base.iterate(dd::DiagDirections{N}, state::Vector{N}) where {N}
    i = 1
    while i <= dd.n && state[i] < 0
        state[i] = -state[i]
        i = i+1
    end
    if i > dd.n
        if dd.n == 1
            # finish here to avoid duplicates
            return nothing
        end
    else
        state[i] = -state[i]
        return (copy(state), state)
    end
end

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

The template can be constructed by passing the dimension. For example, in
dimension two:

```jldoctest dirs_BoxDiag
julia> dirs = BoxDiagDirections(2)
BoxDiagDirections{Float64, Vector{Float64}}(2)

julia> length(dirs) # number of directions
8
```
By default, each direction is represented as a regular vector:

```jldoctest dirs_BoxDiag
julia> eltype(dirs)
Vector{Float64} (alias for Array{Float64, 1})
```
In two dimensions, the directions are normal to the facets of an octagon, i.e.,
the template coincides with [`OctDirections`](@ref).

```jldoctest dirs_BoxDiag
julia> collect(dirs)
8-element Vector{Vector{Float64}}:
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
BoxDiagDirections{Rational{Int64}, Vector{Rational{Int64}}}(10)

julia> length(ans)
1044
```
"""
struct BoxDiagDirections{N, VN} <: AbstractDirections{N, VN}
    n::Int
end

# constructor where only N is specified
BoxDiagDirections{N}(n::Int) where {N} = BoxDiagDirections{N, Vector{N}}(n)

# constructor for type Float64
BoxDiagDirections(n::Int) = BoxDiagDirections{Float64}(n)

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

"""
    PolarDirections{N<:AbstractFloat, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}

Polar directions representation.

### Fields

- `Nφ`    -- length of the partition of the polar angle
- `stack` -- list of computed directions

### Notes

The `PolarDirections` constructor computes a sample of the unit sphere
in ``\\mathbb{R}^2``, which is parameterized by the polar angle
``φ ∈ Dφ := [0, 2π]``; see the Wikipedia entry on the
[polar coordinate system](https://en.wikipedia.org/wiki/Polar_coordinate_system)
for details. The resulting directions are stored in `stack`.

The integer argument ``Nφ`` defines how many samples of ``Dφ`` are taken. The
Cartesian components of each direction are obtained with

```math
[cos(φᵢ), sin(φᵢ)].
```

### Examples

The integer passed as an argument is used to discretize ``φ``:

```jldoctest; filter = r"2246[0-9]*e-16"
julia> pd = PolarDirections(2);

julia> pd.stack
2-element Vector{Vector{Float64}}:
 [1.0, 0.0]
 [-1.0, 1.2246467991473532e-16]

julia> length(pd)
2
```
"""
struct PolarDirections{N<:AbstractFloat, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}
    Nφ::Int
    stack::Vector{VN} # stores the polar directions
end

# constructor where only N is specified
PolarDirections{N}(Nφ::Int) where {N} = PolarDirections{N, Vector{N}}(Nφ)

# constructor for type Float64
PolarDirections(Nφ::Int) = PolarDirections{Float64}(Nφ)

function PolarDirections{N, Vector{N}}(Nφ::Int) where {N}
    if Nφ <= 0
        throw(ArgumentError("Nφ = $Nφ is invalid; it should be at least 1"))
    end
    stack = Vector{Vector{N}}(undef, Nφ)
    φ = range(N(0), stop=N(2*π), length=Nφ+1)  # discretization of the polar angle

    @inbounds for i in 1:Nφ  # skip last (repeated) angle
        stack[i] = N[cos(φ[i]), sin(φ[i])]
    end
    return PolarDirections{N, Vector{N}}(Nφ, stack)
end

Base.eltype(::Type{PolarDirections{N, VN}}) where {N, VN} = VN
Base.length(pd::PolarDirections) = pd.Nφ

# interface functions
dim(pd::PolarDirections) = 2
isbounding(pd::PolarDirections) = pd.Nφ > 2
isnormalized(::Type{<:PolarDirections}) = true

function Base.iterate(pd::PolarDirections{N, Vector{N}}, state::Int=1) where {N}
    state == pd.Nφ + 1 && return nothing
    return (pd.stack[state], state + 1)
end

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
``θ ∈ Dθ := [0, π]`` and ``φ ∈ Dφ := [0, 2π]`` respectively; see the Wikipedia
entry on the
[spherical coordinate system](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
for details.

The integer arguments ``Nθ`` and ``Nφ`` define how many samples along the
domains ``Dθ`` and ``Dφ`` are respectively taken. The Cartesian components of
each direction are obtained with

```math
[sin(θᵢ)*cos(φᵢ), sin(θᵢ)*sin(φᵢ), cos(θᵢ)].
```

The north and south poles are treated separately so that those points
are not considered more than once.

### Examples

The template can be built in different ways. If you pass only one integer, the
same value is used to discretize both ``θ`` and ``φ``:

```jldoctest spherical_directions; filter = r"1232[0-9]*e-17.*2246[0-9]*e-16.*1232[0-9]*e-17"
julia> sd = SphericalDirections(3);

julia> sd.Nθ, sd.Nφ
(3, 3)

julia> length(sd)
4
```

Pass two integers to control the discretization in ``θ`` and in ``φ``
separately:

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

# constructor where only N is specified
SphericalDirections{N}(Nθ::Int, Nφ::Int) where {N} = SphericalDirections{N, Vector{N}}(Nθ::Int, Nφ::Int)

# constructor for type Float64
SphericalDirections(Nθ::Int, Nφ::Int) = SphericalDirections{Float64}(Nθ::Int, Nφ::Int)

# constructor with just one length, interpreted as identical lengths
SphericalDirections(Nθ::Int) = SphericalDirections(Nθ, Nθ)

function SphericalDirections{N, Vector{N}}(Nθ::Int, Nφ::Int) where {N}
    if Nθ <= 1 || Nφ <= 1
        throw(ArgumentError("(Nθ, Nφ) = ($Nθ, $Nφ) is invalid; both should " *
                            "be at least 2"))
    end
    stack = Vector{Vector{N}}()
    θ = range(N(0), stop=N(π), length=Nθ)    # discretization of the azimuthal angle
    φ = range(N(0), stop=N(2*π), length=Nφ)  # discretization of the polar angle

    # add north pole (θ = 0)
    push!(stack, N[0, 0, 1])

    # add south pole (θ = π)
    push!(stack, N[0, 0, -1])

    for φᵢ in φ[1:Nφ-1]  # skip repeated angle
        for θⱼ in θ[2:Nθ-1]  # skip north and south poles
            d = N[sin(θⱼ)*cos(φᵢ), sin(θⱼ)*sin(φᵢ), cos(θⱼ)]
            push!(stack, d)
        end
    end
    return SphericalDirections{N, Vector{N}}(Nθ, Nφ, stack)
end

Base.eltype(::Type{SphericalDirections{N, VN}}) where {N, VN} = VN
Base.length(sd::SphericalDirections) = length(sd.stack)

# interface functions
dim(::SphericalDirections) = 3
isbounding(sd::SphericalDirections) = sd.Nθ > 2 && sd.Nφ > 2
isnormalized(::Type{<:SphericalDirections}) = true

function Base.iterate(sd::SphericalDirections, state::Int=1)
    state == length(sd.stack) + 1 && return nothing
    return (sd.stack[state], state + 1)
end

"""
    CustomDirections{N, VN<:AbstractVector{N}} <: AbstractDirections{N, VN}

User-defined direction vectors.

### Fields

- `directions`          -- list of direction vectors
- `n`                   -- (optional; default: computed from `directions`)
                           dimension
- `check_boundedness`   -- (optional; default: `true`) flag to check boundedness
- `check_normalization` -- (optional; default: `true`) flag to check whether all
                           directions are normalized

### Notes

This struct is a wrapper for a list of user-defined directions. There are fields
for the list of directions, their dimension, and (boolean) cache fields for the
boundedness and normalization properties.
The latter are checked by default upon construction.

To check boundedness, we construct the polyhedron with constraints ``d·x <= 1``
for each direction ``d`` and check if this set is bounded. (Note that the bound
``1`` is arbitrary and that this set may be empty, which however implies
boundedness.)

The dimension will also be determined automatically, unless the empty vector is
passed (in which case the optional argument `n` needs to be specified).

### Examples

Create a template with box directions in dimension two:

```jldoctest
julia> dirs = CustomDirections([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]]);

julia> dirs.directions
4-element Vector{Vector{Float64}}:
 [1.0, 0.0]
 [-1.0, 0.0]
 [0.0, 1.0]
 [0.0, -1.0]

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
    isempty(directions) && throw(ArgumentError("empty direction vectors " *
                                               "need a specified dimension"))
    @inbounds return length(directions[1])
end

function _isbounding(directions::Vector{VN}) where {N, VN<:AbstractVector{N}}
    isempty(directions) && return false

    # check boundedness of the polyhedron `⋂_d d·x <= 1` for directions `d`
    P = HPolyhedron([HalfSpace(dir, one(N)) for dir in directions])
    return isbounded(P)
end

function _isnormalized(directions::Vector{VN}) where {N, VN<:AbstractVector{N}}
    return all(x -> _isapprox(norm(x, 2), one(N)), directions)
end

Base.eltype(::Type{CustomDirections{N, VN}}) where {N, VN} = VN
Base.length(cd::CustomDirections) = length(cd.directions)

# interface functions
dim(cd::CustomDirections) = cd.n
isbounding(cd::Type{<:CustomDirections}) = false  # not a property of the type
isbounding(cd::CustomDirections) = cd.bounded
isnormalized(cd::Type{<:CustomDirections}) = false  # not a property of the type
isnormalized(cd::CustomDirections) = cd.normalized

function Base.iterate(cd::CustomDirections{N}, state::Int=1) where {N}
    if state > length(cd.directions)
        return nothing
    end
    vec = cd.directions[state]
    state = state + 1
    return (vec, state)
end

# instantiation of template for dimension `n` (not possible for all types)
function _get_directions(dir::Type{<:Union{BoxDirections, OctDirections,
                                        DiagDirections, BoxDiagDirections}},
                         n::Int)
    return dir(n)
end

function _get_directions(dir::Type{<:AbstractDirections}, n::Int)
    error("no automatic choice of directions of type $dir possible")
end
