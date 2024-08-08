"""
    VPolygon{N, VN<:AbstractVector{N}} <: AbstractPolygon{N}

Type that represents a polygon by its vertices.

### Fields

- `vertices` -- the list of vertices

### Notes

This type assumes that all vertices are sorted in counter-clockwise fashion.

To ensure this property, the constructor of `VPolygon` runs a convex-hull
algorithm on the vertices by default. This also removes redundant vertices.
If the vertices are known to be sorted, the flag `apply_convex_hull=false` can
be used to skip this preprocessing.

### Examples

A polygon in vertex representation can be constructed by passing the list of
vertices. For example, we can build the right triangle

```jldoctest polygon_vrep
julia> P = VPolygon([[0, 0], [1, 0], [0, 1]]);

julia> P.vertices
3-element Vector{Vector{Int64}}:
 [0, 0]
 [1, 0]
 [0, 1]
```

Alternatively, a `VPolygon` can be constructed passing a matrix of vertices,
where each *column* represents a vertex:

```jldoctest polygon_vrep
julia> M = [0 1 0; 0 0 1.]
2×3 Matrix{Float64}:
 0.0  1.0  0.0
 0.0  0.0  1.0

julia> P = VPolygon(M);

julia> P.vertices
3-element Vector{Vector{Float64}}:
 [0.0, 0.0]
 [1.0, 0.0]
 [0.0, 1.0]
```
"""
struct VPolygon{N,VN<:AbstractVector{N}} <: AbstractPolygon{N}
    vertices::Vector{VN}

    # default constructor that applies a convex hull algorithm
    function VPolygon(vertices::Vector{VN};
                      apply_convex_hull::Bool=true,
                      algorithm::String="monotone_chain") where {N,VN<:AbstractVector{N}}
        if apply_convex_hull
            vertices = convex_hull(vertices; algorithm=algorithm)
        end
        return new{N,VN}(vertices)
    end
end

# constructor with empty vertices list
VPolygon{N}() where {N} = VPolygon(Vector{Vector{N}}(); apply_convex_hull=false)

# constructor with no vertices of type Float64
VPolygon() = VPolygon{Float64}()

# constructor from rectangular matrix
function VPolygon(vertices_matrix::MT; apply_convex_hull::Bool=true,
                  algorithm::String="monotone_chain") where {MT<:AbstractMatrix}
    @assert size(vertices_matrix, 1) == 2 "the number of rows of the matrix " *
                                          "of vertices should be 2, but it is $(size(vertices_matrix, 1))"

    vertices = _to_colVector(vertices_matrix)
    return VPolygon(vertices; apply_convex_hull=apply_convex_hull,
                    algorithm=algorithm)
end
