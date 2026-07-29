module GeometryBasicsExt

using LazySets: LazySet, polyhedron, triangulate_faces, _area_triangle_3D!,
                @validate
using GeometryBasics: coordinates, faces
using Polyhedra: Mesh
import LazySets: triangulate_faces, _area_polytope_3D

function _area_polytope_3D(P::LazySet)
    points, connections = triangulate_faces(P)
    N = (eltype(P) <: AbstractFloat) ? eltype(P) : Float64  # `_area_triangle_3D!` uses `sqrt`
    res = zero(N)
    M = ones(N, 3, 3)
    @inbounds for triple in connections
        x, y, z = points[:, triple[1]], points[:, triple[2]], points[:, triple[3]]
        res += _area_triangle_3D!(M, x, y, z)
    end
    return res
end

"""
    triangulate_faces(P::LazySet)

Triangulate the faces of a three-dimensional polytopic set.

### Input

- `P` -- three-dimensional polytopic set

### Output

A tuple `(p, c)` where `p` is a matrix, with each column containing a point, and
`c` is a list of 3-tuples containing the indices of the points in each triangle.

### Notes

This function triangulates all faces of a 3D polytope. The result is a list of (flat)
triangles in 3D which describe the boundary of `P`.

`P` must contain at least three vertices.
"""
@validate function triangulate_faces(P::LazySet)
    mes = Mesh(polyhedron(P))
    coords = coordinates(mes)
    connection = faces(mes)

    ntriangles = length(connection)
    npoints = length(coords)
    @assert npoints == 3 * ntriangles "each triangle should have 3 vertices"
    points = Matrix{Float32}(undef, 3, npoints)

    for i in 1:npoints
        points[:, i] .= coords[i].data
    end

    connection_tup = getfield.(connection, :data)

    return points, connection_tup
end

end  # module
