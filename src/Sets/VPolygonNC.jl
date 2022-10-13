export VPolygonNC

"""
    VPolygonNC{N, VN<:AbstractVector{N}} <: LazySet{N}

Type that represents a non-convex polygon by its vertices.

### Fields

- `vertices` -- the list of vertices (in clockwise or counter-clockwise order)

### Examples

A non-convex polygon in vertex representation can be constructed by passing the
list of vertices. For example, we can build a tooth:

```jldoctest polygon_v_ncrep
julia> P = VPolygonNC([[0.0, 0], [0, 2], [2, 2], [2, 0], [1, 1]]);

julia> P.vertices
5-element Vector{Vector{Float64}}:
 [0.0, 0.0]
 [0.0, 2.0]
 [2.0, 2.0]
 [2.0, 0.0]
 [1.0, 1.0]
```
"""
struct VPolygonNC{N, VN<:AbstractVector{N}} <: LazySet{N}
    vertices::Vector{VN}

    function VPolygonNC(vertices::Vector{VN}) where {N, VN<:AbstractVector{N}}
        return new{N, VN}(vertices)
    end
end

# constructor with empty vertices list
VPolygonNC{N}() where {N} = VPolygonNC(Vector{Vector{N}}())

# constructor with no vertices of type Float64
VPolygonNC() = VPolygonNC{Float64}()

function dim(P::VPolygonNC)
    return 2
end

function isbounded(P::VPolygonNC)
    return true
end

function vertices_list(P::VPolygonNC)
    return P.vertices
end

function plot_recipe(P::VPolygonNC{N}, ε=zero(N)) where {N}
    vlist = vertices_list(P)
    return _plot_recipe_2d_vlist(vlist, N)
end
