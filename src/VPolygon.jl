import Base.<=

# Polygons in vertex representation
export VPolygon, vertices_list, singleton_list

"""
    VPolygon <: LazySet

Type that represents a polygon by its vertices.

### Fields

- `vertices_list` -- the list of vertices

### Notes

The constructor of `VPolygon` runs a convex hull algorithm, and the given vertices
are sorted in counter-clockwise fashion. If you don't want to take the
convex hull, set the `apply_convex_hull=false` flag when instantiating the constructor.
"""
struct VPolygon <: LazySet
    vertices_list::Vector{Vector{Float64}}

    function VPolygon(vertices_list; apply_convex_hull=true, algorithm="monotone_chain")
        if apply_convex_hull
            return new(convex_hull(vertices_list, algorithm=algorithm))
        else
            return new(vertices_list)
        end
    end
end
VPolygon() = VPolygon([])

"""
    dim(P)

Return the ambient dimension of the polygon.

### Input

- `P` -- polygon in vertex representation
"""
function dim(P::VPolygon)
    2
end

"""
    vertices_list(P)

Return the list of vertices of a convex polygon in vertex representation.

### Input

- `P` -- a polygon given in vertex representation

### Output

List of vertices as an array of vertex pairs, `Vector{Vector{Float64}}`.
"""
function vertices_list(P::VPolygon)::Vector{Vector{Float64}}
    return P.vertices_list
end

"""
    singleton_list(P)

Return the vertices of a convex polygon in vertex representation as a list of
singletons.

### Input

- `P` -- a polygon given in vertex representation

### Output

List of vertices as an array of vertex pairs, `Vector{Singleton{Float64}}`.
"""
function singleton_list(P::VPolygon)
    return [Singleton(vi) for vi in P.vertices_list]
end
