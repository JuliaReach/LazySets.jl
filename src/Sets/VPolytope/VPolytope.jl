"""
    VPolytope{N, VN<:AbstractVector{N}, VT<:AbstractVector{VN}} <: AbstractPolytope{N}

Type that represents a convex polytope in vertex representation.

### Fields

- `vertices` -- list of vertices

### Examples

A polytope in vertex representation can be constructed by passing the list of
vertices. For example, we can build the tetrahedron:

```jldoctest
julia> P = VPolytope([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]);

julia> P.vertices
4-element Vector{Vector{Int64}}:
 [0, 0, 0]
 [1, 0, 0]
 [0, 1, 0]
 [0, 0, 1]
```

Alternatively, a `VPolytope` can be constructed passing a matrix of vertices,
where each *column* represents a vertex:

```jldoctest
julia> M = [0 0 0; 1 0 0; 0 1 0; 0 0 1]'
3×4 adjoint(::Matrix{Int64}) with eltype Int64:
 0  1  0  0
 0  0  1  0
 0  0  0  1

julia> P = VPolytope(M);

julia> P.vertices
4-element Vector{Vector{Int64}}:
 [0, 0, 0]
 [1, 0, 0]
 [0, 1, 0]
 [0, 0, 1]
```
"""
struct VPolytope{N,VN<:AbstractVector{N},VT<:AbstractVector{VN}} <: AbstractPolytope{N}
    vertices::VT
end

# constructor with empty vertices list
VPolytope{N}() where {N} = VPolytope(Vector{Vector{N}}())

# constructor with no vertices of type Float64
VPolytope() = VPolytope{Float64}()

# constructor from rectangular matrix
function VPolytope(vertices_matrix::MT) where {N,MT<:AbstractMatrix{N}}
    vertices = [vertices_matrix[:, j] for j in axes(vertices_matrix, 2)]
    return VPolytope(vertices)
end

function load_Polyhedra_VPolytope()
    return quote
        using .Polyhedra: VRep

        # VPolytope from a VRep
        function VPolytope(P::VRep)
            vertices = collect(Polyhedra.points(P))
            return VPolytope(vertices)
        end
    end
end  # load_Polyhedra_VPolytope
