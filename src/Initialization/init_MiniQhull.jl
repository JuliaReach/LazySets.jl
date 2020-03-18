using .MiniQhull
import .MiniQhull: delaunay
export delaunay

"""
    delaunay(X::LazySet)

Compute the Delaunay triangulation of the given convex set.

### Input

- `X` -- set

### Output

A union of polytopes in vertex representation.

### Notes

This function requires that you have properly installed
the Julia package [MiniQhull.jl](https://github.com/gridap/MiniQhull.jl), including
the library [Qhull](http://www.qhull.org/).

The method works in arbitrary dimension and the requirement is that the list of
vertices of `X` can be obtained.
"""
function delaunay(X::LazySet)
    n = dim(X)
    v = vertices_list(X)
    m = length(v)
    coordinates = vcat(v...)
    connectivity_matrix = delaunay(n, m, coordinates)
    nelements = size(connectivity_matrix, 2)
    elements = [VPolytope(v[connectivity_matrix[:, i]]) for i in 1:nelements]
    return UnionSetArray(elements)
end
