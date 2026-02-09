# The algorithm adds an edge for each consecutive pair of vertices.
# Since the vertices are already ordered in counter-clockwise fashion (CCW), the
# constraints will be sorted (CCW) as well.
function constraints_list(P::VPolygon)
    vl = P.vertices
    m = length(vl)
    if m == 0
        # no vertex
        N = eltype(P)
        clist = _infeasible_constraints_list(2; N=N)
    elseif m == 1
        # only one vertex -> use function for singletons
        require(@__MODULE__, :LazySets; fun_name="convert")
        clist = _constraints_list_singleton_Vector(vl[1])
    elseif m == 2
        # only two vertices -> use function for line segments
        require(@__MODULE__, :LazySets; fun_name="convert")
        clist = constraints_list(LineSegment(vl[1], vl[2]))
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
    end
    return clist
end
