import Base.<=

# Polygons in vertex representation
export VPolygon, tovrep, vertices_list, is_contained

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

    function VPolygon(vertices_list; apply_convex_hull=true, algorithm="andrew_monotone_chain")
        if apply_convex_hull
            return new(convex_hull(vertices_list, algorithm=algorithm))
        else
            return new(vertices_list)
        end
    end
end
VPolygon() = VPolygon([])


"""
    vertices_list(P)

Return the list of vertices of a convex polygon in vertex representation.

### Input

- `P` -- a polygon given in vertex representation

### Output

List of vertices as an array of vertex pairs, Vector{Vector{Float64}}.
"""
function vertices_list(P::VPolygon)::Vector{Vector{Float64}}
    return P.vertices_list
end

"""
    plot_polygon(P::VPolygon; ...)

Plot a polygon given in vertex representation.

### Input

- `P` -- a polygon in vertex representation

### Examples

```julia
julia> using LazySets, Plots
julia> P = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]])
julia> plot(P)
```
"""
@recipe function plot_polygon(P::VPolygon;
                              color="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]

     x, y
end

"""
    plot_polygons(P::Vector{VPolygon}; ...)

Plot an array of polygons given in vertex representation.

### Input

- `P` -- an array of polygons in vertex representation

### Examples

```julia
julia> using LazySets, Plots
julia> P1 = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]])
julia> P2 = VPolygon([[0.3, 0.3], [0.2, 0.3], [0.2, 0.2], [0.3, 0.2]])
julia> plot([P1, P2])
```
"""
@recipe function plot_polygons(P::Vector{VPolygon};
                              seriescolor="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    for Pi in P
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end
