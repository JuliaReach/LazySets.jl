"""
    tohrep(P::VPolygon, ::Type{HPOLYGON}=HPolygon) where {HPOLYGON<:AbstractHPolygon}

Build a constraint representation of the given polygon.

### Input

- `P`        -- polygon in vertex representation
- `HPOLYGON` -- (optional, default: `HPolygon`) type of target polygon

### Output

A polygon in constraint representation, an `AbstractHPolygon`.

### Algorithm

The algorithm adds an edge for each consecutive pair of vertices.
Since the vertices are already ordered in counter-clockwise fashion (CCW), the
constraints will be sorted automatically (CCW).
"""
function tohrep(P::VPolygon, ::Type{HPOLYGON}=HPolygon) where {HPOLYGON<:AbstractHPolygon}
    vl = P.vertices
    n = length(vl)
    if n == 0
        # no vertex
        N = eltype(P)
        constraints_list = _infeasible_constraints_list(2; N=N)
    elseif n == 1
        # only one vertex -> use function for singletons
        require(@__MODULE__, :LazySets; fun_name="convert")
        return convert(HPOLYGON, Singleton(vl[1]))
    elseif n == 2
        # only two vertices -> use function for line segments
        require(@__MODULE__, :LazySets; fun_name="convert")
        return convert(HPOLYGON, LineSegment(vl[1], vl[2]))
    else
        # find right-most vertex
        i = div(n, 2)
        x = vl[i][1]
        while i > 1 && vl[i - 1][1] > x
            # search forward in list
            i = i - 1
            x = vl[i][1]
        end
        while i < n && vl[i + 1][1] > x
            # search backward in list
            i = i + 1
            x = vl[i][1]
        end

        # create constraints ordered in CCW starting at the right-most index
        upper_hull = [halfspace_left(vl[j], vl[j + 1]) for j in i:(length(vl) - 1)]
        mid_hull = [halfspace_left(vl[end], vl[1])]
        lower_hull = [halfspace_left(vl[j], vl[j + 1]) for j in 1:(i - 1)]
        constraints_list = vcat(upper_hull, mid_hull, lower_hull)
    end
    return HPOLYGON(constraints_list)
end
