export Polygon

"""
    Polygon{N, VN<:AbstractVector{N}} <: LazySet{N}

Type that represents a non-convex polygon by its vertices.

### Fields

- `vertices` -- the list of vertices (in clockwise or counter-clockwise order)

### Examples

A non-convex polygon in vertex representation can be constructed by passing the
list of vertices. For example, we can build a tooth:

```jldoctest polygon_v_ncrep
julia> P = Polygon([[0.0, 0], [0, 2], [2, 2], [2, 0], [1, 1]]);

julia> P.vertices
5-element Vector{Vector{Float64}}:
 [0.0, 0.0]
 [0.0, 2.0]
 [2.0, 2.0]
 [2.0, 0.0]
 [1.0, 1.0]
```
"""
struct Polygon{N,VN<:AbstractVector{N}} <: LazySet{N}
    vertices::Vector{VN}

    function Polygon(vertices::Vector{VN}) where {N,VN<:AbstractVector{N}}
        return new{N,VN}(vertices)
    end
end

# constructor with empty vertices list
Polygon{N}() where {N} = Polygon(Vector{Vector{N}}())

# constructor with no vertices of type Float64
Polygon() = Polygon{Float64}()

function isoperationtype(P::Type{<:Polygon})
    return false
end

function isconvextype(P::Type{<:Polygon})
    return false
end

function isboundedtype(P::Type{<:Polygon})
    return true
end

function dim(P::Polygon)
    return 2
end

function isbounded(P::Polygon)
    return true
end

function plot_recipe(P::Polygon{N}, ε=zero(N)) where {N}
    vlist = P.vertices
    return _plot_recipe_2d_vlist(vlist, N)
end
