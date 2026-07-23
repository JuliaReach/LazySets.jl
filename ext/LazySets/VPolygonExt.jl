using Base: extrema
using LazySets: vertices_list, permute, convex_hull!, halfspace_left,
                _constraints_list_singleton_Vector,
                _infeasible_constraints_list, _intersection_vrep_2d,
                _sort_constraints, @validate
using LazySets.EmptySetModule: EmptySet
using LazySets.IntervalModule: Interval
using LazySets.VPolygonModule: VPolygon
using LazySets.LineSegmentModule: LineSegment
import LazySets.API: constraints_list, project, convex_hull, intersection
import LazySets: remove_redundant_vertices!

# The algorithm adds an edge for each consecutive pair of vertices.
# Since the vertices are already ordered in counter-clockwise fashion (CCW), the
# constraints will be sorted (CCW) as well.
function constraints_list(P::VPolygon; sort_constraints::Bool=false)
    vl = P.vertices
    m = length(vl)
    if m == 0
        # no vertex
        N = eltype(P)
        clist = _infeasible_constraints_list(2; N=N)
    elseif m == 1
        # only one vertex -> use function for singletons
        clist = _constraints_list_singleton_Vector(vl[1])
    elseif m == 2
        # only two vertices -> use function for line segments
        clist = constraints_list(LineSegment(vl[1], vl[2]); sort_constraints)
    else
        # find right-most vertex
        i = div(m, 2)
        x = vl[i][1]
        while i > 1 && vl[i - 1][1] > x
            # search forward in list
            i = i - 1
            x = vl[i][1]
        end
        while i < m && vl[i + 1][1] > x
            # search backward in list
            i = i + 1
            x = vl[i][1]
        end

        # create constraints ordered in CCW starting at the right-most index
        upper_hull = [halfspace_left(vl[j], vl[j + 1]) for j in i:(length(vl) - 1)]
        mid_hull = [halfspace_left(vl[end], vl[1])]
        lower_hull = [halfspace_left(vl[j], vl[j + 1]) for j in 1:(i - 1)]
        clist = vcat(upper_hull, mid_hull, lower_hull)

        if sort_constraints
            # result is not sorted according to `⪯` yet, so sort it
            # TODO can the code be changed to get rid of the sorting? here is a counterexample:
            # `VPolygon{Float64, Vector{Float64}}([[0.0, 0.0], [-1.0, -1.0], [1.0, -1.0]])`
            # julia> constraints_list(R; sort_constraints=true)
            # 3-element Vector{HalfSpace{Float64, Vector{Float64}}}:
            # HalfSpace{Float64, Vector{Float64}}([1.0, 1.0], 0.0)
            # HalfSpace{Float64, Vector{Float64}}([-1.0, 1.0], 0.0)
            # HalfSpace{Float64, Vector{Float64}}([0.0, -2.0], 2.0)
            # julia> constraints_list(R; sort_constraints=false)
            # 3-element Vector{HalfSpace{Float64, Vector{Float64}}}:
            # HalfSpace{Float64, Vector{Float64}}([-1.0, 1.0], 0.0)
            # HalfSpace{Float64, Vector{Float64}}([0.0, -2.0], 2.0)
            # HalfSpace{Float64, Vector{Float64}}([1.0, 1.0], 0.0)
            clist = _sort_constraints(clist)
        end
    end
    return clist
end

@validate function project(V::VPolygon, block::AbstractVector{Int}; kwargs...)
    if length(block) == 1
        l, h = extrema(V, block[1])
        return Interval(l, h)
    end
    # length(block) == 2
    if block[1] == 1 && block[2] == 2
        return V  # no projection
    else
        # block[1] == 2 && block[2] == 1
        return permute(V, block)  # swap dimensions
    end
end

# TODO outsource to helper function in LazySets
@validate function convex_hull(P::VPolygon, Q::VPolygon;
                               algorithm::String="monotone_chain")
    vunion = [P.vertices; Q.vertices]
    convex_hull!(vunion; algorithm=algorithm)
    return VPolygon(vunion; apply_convex_hull=false)
end

"""
# Extended help

    intersection(P1::VPolygon, P2::VPolygon; apply_convex_hull::Bool=true)

### Output

A `VPolygon`, or an `EmptySet` if the intersection is empty.

### Algorithm

This function applies the [Sutherland–Hodgman polygon clipping
algorithm](https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm).
The implementation is based on the one found in
[rosetta code](http://www.rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Julia).
"""
@validate function intersection(P1::VPolygon, P2::VPolygon; apply_convex_hull::Bool=true)
    v1 = vertices_list(P1)
    v2 = vertices_list(P2)
    v12 = _intersection_vrep_2d(v1, v2)

    if isempty(v12)
        N = promote_type(eltype(P1), eltype(P2))
        return EmptySet{N}(2)
    else
        return VPolygon(v12; apply_convex_hull=apply_convex_hull)
    end
end

# TODO move back to VPolygonModule once `convex_hull!` is an API function
"""
    remove_redundant_vertices!(P::VPolygon;
                               [algorithm]::String="monotone_chain")

Remove the redundant vertices from the given polygon in-place.

### Input

- `P`         -- polygon in vertex representation
- `algorithm` -- (optional, default: "monotone_chain") the algorithm used to
                 compute the convex hull

### Output

The modified polygon whose redundant vertices have been removed.

### Algorithm

A convex-hull algorithm is used to compute the convex hull of the vertices of
the polygon `P`; see `?convex_hull` for details on the available algorithms.
The vertices are sorted in counter-clockwise fashion.
"""
function remove_redundant_vertices!(P::VPolygon;
                                    algorithm::String="monotone_chain")
    convex_hull!(P.vertices; algorithm=algorithm)
    return P
end
