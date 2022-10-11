export VPolygonNC

"""
    VPolygonNC{N, VN<:AbstractVector{N}} <: LazySet{N}

Type that represents a non-convex polygon by its vertices.

### Fields

- `vertices` -- the list of vertices

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
