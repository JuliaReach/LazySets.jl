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

# constructor from rectangular matrix
function Tetrahedron(vertices_matrix::MT) where {N,MT<:AbstractMatrix{N}}
    vertices = collect(eachcol(vertices_matrix))
    return Tetrahedron(vertices)
end
