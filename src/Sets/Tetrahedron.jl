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

julia> T = Tetrahedron(vertices);

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

    function Tetrahedron(vertices::VT) where {N,VN<:AbstractVector{N},VT<:AbstractVector{VN}}
        @assert length(vertices) == 4 "a tetrahedron requires four vertices"
        return new{N,VN,VT}(vertices)
    end
end

isoperationtype(::Type{<:Tetrahedron}) = false

# constructor with empty vertices list
Tetrahedron{N}() where {N} = Tetrahedron(Vector{Vector{N}}())

# constructor with no vertices of type Float64
Tetrahedron() = Tetrahedron{Float64}()

# constructor from rectangular matrix
function Tetrahedron(vertices_matrix::MT) where {N,MT<:AbstractMatrix{N}}
    vertices = collect(eachcol(vertices_matrix))
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

- `x` -- point/vector
- `P` -- tetrahedron in vertex representation

### Output

`true` iff ``x ∈ P``.

### Algorithm

For each plane of the tetrahedron, we check if the point `x` is on the same side as the remaining vertex.
We need to check this for each plane [1].

[1] https://stackoverflow.com/questions/25179693/how-to-check-whether-the-point-is-in-the-tetrahedron-or-not
"""
function ∈(x::AbstractVector, T::Tetrahedron)
    v = T.vertices
    return same_side(v[1], v[2], v[3], v[4], x) &&
           same_side(v[4], v[1], v[2], v[3], x) &&
           same_side(v[3], v[4], v[1], v[2], x) &&
           same_side(v[2], v[3], v[4], v[1], x)
end

# Return `true` iff point `x` lies in the same half-space as `v4` with respect to the hyperplane determined by `v1`, `v2` and `v3`.
function same_side(v1, v2, v3, v4, x)
    normal = cross(v2 - v1, v3 - v1)
    dotv4 = dot(normal, v4 - v1)
    dotx = dot(normal, x - v1)
    return signbit(dotv4) == signbit(dotx)
end

function rand(::Type{Tetrahedron}; N::Type{<:Real}=Float64, rng::AbstractRNG=GLOBAL_RNG, seed::Union{Int,Nothing}=nothing)
    P = rand(VPolytope; N=N, dim=3, rng=rng, seed=seed, num_vertices=4)
    vertices = P.vertices
    return Tetrahedron(vertices)
end

function constraints_list(T::Tetrahedron)
    constraints_list(convert(VPolytope, T))
end
