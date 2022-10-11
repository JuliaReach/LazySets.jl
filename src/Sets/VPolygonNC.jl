export VPolygonNC

"""
    VPolygonNC{N, VN<:AbstractVector{N}} <: LazySet{N}

Type that represents a non-convex polygon by its vertices.

### Fields

- `vertices` -- the list of vertices

### Examples

A non-convex polygon in vertex representation can be constructed by passing the list of
vertices. For example, we can build a teeth:

```jldoctest polygon_v_ncrep
julia> NCP = VPolygonNC( [[0., 0], [0, 2], [2, 2], [2, 0], [1, 1]]);

julia> NCP.vertices
5-element Vector{Vector{Float64}}:
 [0.0, 0.0]
 [0.0, 2.0]
 [2.0, 2.0]
 [2.0, 0.0]
 [1.0, 1.0]
```

```
"""
struct VPolygonNC{N, VN<:AbstractVector{N}} <: LazySet{N}
    vertices::Vector{VN}

    # default constructor that applies a convex hull algorithm
    function VPolygonNC(vertices::Vector{VN}) where {N, VN<:AbstractVector{N}}
        return new{N, VN}(vertices)
    end
end

# constructor with empty vertices list
VPolygonNC{N}() where {N} = VPolygonNC(Vector{Vector{N}}())

# constructor with no vertices of type Float64
VPolygonNC() = VPolygonNC{Float64}()
