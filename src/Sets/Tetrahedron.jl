export Tetrahedron

"""
    Tetrahedron{N, VN<:AbstractVector{N}, VT<:AbstractVector{VN}} <: AbstractPolytope{N}

Type that represents a tetrahedron in vertex representation.

### Fields

- `vertices` -- list of vertices

### Examples

A tetrahedron can be constructed by passing the list of vertices.
The following builds the tetrahedron with edge length 2 from the [wikipedia page Tetrahedron](https://en.wikipedia.org/wiki/Tetrahedron): 

```jldoctest
julia> vertices = [[1, 0, -1/sqrt(2)], [-1, 0, -1/sqrt(2)], [0, 1, 1/sqrt(2)], [0, -1, 1/sqrt(2)]];

julia> T = Tetrahedron(vertices); T isa Tetrahedron
true

julia> dim(T)
3

julia> zeros(3) ∈ T
true

julia> σ(ones(3), T)
3-element Vector{Float64}:
 0.0
 1.0
 0.7071067811865475
```
"""
struct Tetrahedron{N,VN<:AbstractVector{N},VT<:AbstractVector{VN}} <: AbstractPolytope{N}
    vertices::VT
end

isoperationtype(::Type{<:Tetrahedron}) = false

# constructor with empty vertices list
Tetrahedron{N}() where {N} = Tetrahedron(Vector{Vector{N}}())

# constructor with no vertices of type Float64
Tetrahedron() = Tetrahedron{Float64}()

# constructor from rectangular matrix
function Tetrahedron(vertices_matrix::MT) where {N,MT<:AbstractMatrix{N}}
    vertices = [vertices_matrix[:, j] for j in axes(vertices_matrix, 2)]
    return Tetrahedron(vertices)
end

function dim(::Tetrahedron)
    return 3
end

"""
    σ(d::AbstractVector, P::VPolytope)

Return a support vector of a tetrahedron in a given direction.

### Input

- `d` -- direction
- `P` -- tetrahedron

### Output

A support vector in the given direction.

### Algorithm

Currently falls back to the `VPolytope` implementation.
"""
function σ(d::AbstractVector, T::Tetrahedron)
    return σ(d, convert(VPolytope, T))
end

"""
    ∈(x::AbstractVector, T::Tetrahedron)

Check whether a given point is contained in a tetrahedron.

### Input

- `x`      -- point/vector
- `P`      -- tetrahedron in vertex representation

### Output

`true` iff ``x ∈ P``.

### Algorithm

Currently falls back to the `VPolytope` implementation.
TODO Implement efficient version, see LazySets#3259.
"""
function ∈(x::AbstractVector, T::Tetrahedron)
    return x ∈ convert(VPolytope, T)
end

function rand(::Type{Tetrahedron}; N::Type{<:Real}=Float64, rng::AbstractRNG=GLOBAL_RNG, seed::Union{Int, Nothing}=nothing)
    # Dummy result for tests.
    vertices = [N[1, 0, -1/sqrt(2)], N[-1, 0, -1/sqrt(2)], N[0, 1, 1/sqrt(2)], N[0, -1, 1/sqrt(2)]];
    return Tetrahedron(vertices)
end

function constraints_list(T::Tetrahedron)
    constraints_list(convert(VPolytope, T))
end
