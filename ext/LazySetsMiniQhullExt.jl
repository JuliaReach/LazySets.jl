module LazySetsMiniQhullExt

using MiniQhull: delaunay
import LazySets: EmptySet, LazySet, UnionSetArray, VPolytope, dim, ispolytopic,
                 plot_recipe, vertices_list, _plot_recipe_3d_polytope,
                 _triangulate_delaunay

function _triangulate_delaunay(X::LazySet; compute_triangles_3d::Bool=false)
    vlist, connect_mat = _delaunay_vlist_connectivity(X; compute_triangles_3d=compute_triangles_3d)
    nsimplices = size(connect_mat, 2)
    if compute_triangles_3d
        simplices = [VPolytope(vlist[connect_mat[1:3, j]]) for j in 1:nsimplices]
    else
        simplices = [VPolytope(vlist[connect_mat[:, j]]) for j in 1:nsimplices]
    end
    return UnionSetArray(simplices)
end

# compute the vertices and the connectivity matrix of the Delaunay triangulation
#
# if the flag `compute_triangles_3d` is set, the resulting matrix still has four
# rows, but the last row has no meaning
function _delaunay_vlist_connectivity(X::LazySet; compute_triangles_3d::Bool=false)
    n = dim(X)
    @assert !compute_triangles_3d || n == 3 "the `compute_triangles_3d` " *
                                            "option requires 3D inputs"
    vlist = vertices_list(X)
    m = length(vlist)
    coordinates = vcat(vlist...)
    flags = compute_triangles_3d ? "qhull Qt" : nothing
    connectivity_matrix = delaunay(n, m, coordinates, flags)
    return vlist, connectivity_matrix
end

function _plot_recipe_3d_polytope(P::LazySet, N=eltype(P))
    @assert ispolytopic(P) "3D plotting is only available for polytopes"

    vlist, C = _delaunay_vlist_connectivity(P; compute_triangles_3d=true)

    m = length(vlist)
    if m == 0
        @warn "received a polyhedron with no vertices during plotting"
        return plot_recipe(EmptySet{N}(2), zero(N))
    end

    x = Vector{N}(undef, m)
    y = Vector{N}(undef, m)
    z = Vector{N}(undef, m)
    @inbounds for (i, vi) in enumerate(vlist)
        x[i] = vi[1]
        y[i] = vi[2]
        z[i] = vi[3]
    end

    l = size(C, 2)
    i = Vector{Int}(undef, l)
    j = Vector{Int}(undef, l)
    k = Vector{Int}(undef, l)
    @inbounds for idx in 1:l
        # normalization: -1 for zero indexing; convert to Int on 64-bit systems
        i[idx] = Int(C[1, idx] - 1)
        j[idx] = Int(C[2, idx] - 1)
        k[idx] = Int(C[3, idx] - 1)
    end

    return x, y, z, i, j, k
end

end  # module
