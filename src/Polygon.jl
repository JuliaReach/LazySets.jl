"""
    HPolygon <: LazySet

Type that represents a convex polygon (in H-representation).

FIELDS:

- ``constraints`` --  an array of linear constraints

EXMPLES:

A triangle in the first quadrant::

    julia> HPolygon()
"""
mutable struct HPolygon <: LazySet
    constraints::Array{LinearConstraint, 1}
end
HPolygon() = HPolygon([])

function dim(P::HPolygon)::Int64
    2
end

import Base.<=

"""
    u <= v

States if arg(u) [2π] <= arg(v) [2π].

INPUT:

- ``u`` --  a first direction
- ``v`` --  a second direction

OUTPUT:

True iff arg(u) [2π] <= arg(v) [2π]

NOTE:

The argument is measured in counter-clockwise fashion, with the 0 being the
direction (1, 0).
"""
function <=(u::Union{Vector{Float64}, SparseVector{Float64,Int64}},
            v::Union{Vector{Float64}, SparseVector{Float64,Int64}})
    return jump2pi(atan2(u[2], u[1])) <= jump2pi(atan2(v[2], v[1]))
end

"""
    addconstraint!(p, c)

Add a linear constraint to a polygon keeping the constraints sorted by their
normal directions.

INPUT:

- ``p`` -- a polygon
- ``c`` -- the linear constraint to add
"""
function addconstraint!(p::HPolygon, c::LinearConstraint)
    i = length(p.constraints)
    while i > 0 && c.a <= p.constraints[i].a
        i -= 1
    end
    # here p.constraints[i] < c
    insert!(p.constraints, i+1, c)
end

"""
    σ(d, P)

Return the support vector of a polygon in a given direction.
Return the zero vector if there are no constraints (i.e., the universal
polytope).

INPUT:

- ``d`` -- direction
- ``p`` -- polyhedron in H-representation
"""
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, p::HPolygon)::Vector{Float64}
    n = length(p.constraints)
    if n == 0
        error("this polygon is empty")
    end
    i = 1
    while i <= n && p.constraints[i].a <= d
        i += 1
    end
    if i == 1 || i == n+1
        intersection(Line(p.constraints[1]), Line(p.constraints[n]))
    else
        intersection(Line(p.constraints[i]), Line(p.constraints[i-1]))
    end
end


"""
    is_contained(x, P)

Return whether a given vector is contained in the polygon.

INPUT:

- ``x`` -- vector
- ``P`` -- polygon

OUTPUT:

::Bool : true iff x ∈ P
"""
function is_contained(x::Vector{Float64}, P::HPolygon)::Bool
    if (length(x) != 2)
        false
    else
        res = true
        for c in P.constraints
            res = res && dot(c.a, x) <= c.b
        end
        res
    end
end

"""
    HPolygonOpt <: LazySet

Type that represents a convex polygon (in H-representation).
This structure is optimized to evaluate the support function/vector with a large
sequence of directions, which are one to one close.

FIELDS:

- ``P``   -- polygon
- ``ind`` -- an index in the list of constraints to begin the search to
             evaluate the support functions.
"""
mutable struct HPolygonOpt <: LazySet
    constraints::Array{LinearConstraint, 1}
    ind::Int64
    HPolygonOpt(constraints) = new(constraints, 1)
    HPolygonOpt(constraints, ind) = new(constraints, ind)
end
HPolygonOpt(H::HPolygon) = HPolygonOpt(H.constraints)

"""
    dim(P)

Return the ambient dimension of the optimized polygon.

INPUT :

- ``P`` -- optimized polyhedron in H-representation
"""
function dim(P::HPolygonOpt)::Int64
    return 2
end

"""
    σ(d, P)

Return the support vector of the optimized polygon in a given direction.

INPUT:

- ``d`` -- direction
- ``P`` -- polyhedron in H-representation
"""
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, p::HPolygonOpt)::Vector{Float64}
    n = length(p.constraints)
    cl = p.constraints
    if (d <= cl[p.ind].a)
        k = p.ind-1
        while (k >= 1 && d <= cl[k].a)
            k -= 1
        end
        if (k == 0)
            p.ind = n
            return intersection(Line(cl[n]), Line(cl[1]))
        else
            p.ind = k
            return intersection(Line(cl[k]), Line(cl[k+1]))
        end
    else
        k = p.ind+1
        while (k <= n && cl[k].a <= d)
            k += 1
        end
        if (k == n+1)
            p.ind = n
            return intersection(Line(cl[n]), Line(cl[1]))
        else
            p.ind = k-1
            return intersection(Line(cl[k-1]), Line(cl[k]))
        end
    end
end

"""
    is_contained(x, P)

Return whether a given vector is contained in the optimized polygon.

INPUT:

- ``x`` -- vector
- ``P`` -- polygon

OUTPUT:

True iff x ∈ P
"""
function is_contained(x::Vector{Float64}, P::HPolygonOpt)::Bool
    if (length(x) != 2)
        false
    else
        res = true
        for c in P.constraints
            res = res && dot(c.a, x) <= c.b
        end
        res
    end
end


"""
    VPolygon

Type that represents a polygon by its vertices.

FIELDS:

- ``vl`` -- the list of vertices
"""
mutable struct VPolygon
    vl::Array{Vector{Float64}, 1}
end
VPolygon() = VPolygon([])

"""
    tovrep(s)

Build a vertex representation of the given polygon.

INPUT:

- ``s`` -- a polygon in H-representation, HPolygon. The linear constraints are
           assumed sorted by their normal directions.

OUTPUT:

The same polygon in a vertex representation, VPolygon.
"""
function tovrep(s::HPolygon)
    n = length(s.constraints)
    p = VPolygon()
    for i in 1:n-1
        cur = intersection(Line(s.constraints[i]), Line(s.constraints[i+1]))
        push!(p.vl, cur)
    end
    cur = intersection(Line(s.constraints[n]), Line(s.constraints[1]))
    push!(p.vl, cur)
    p
end

"""
    tovrep(po)

Build a vertex representation of the given polygon.

INPUT :

- ``po`` -- a polygon in H-representation. The linear constraints are
            assumed sorted by their normal directions.

OUTPUT :

The same polygon in a vertex representation.
"""
function tovrep(po::HPolygonOpt)
    tovrep(po.p)
end

"""
    vertices_list(po::Union{HPolygon, HPolygonOpt})

Return the list of vertices of a convex polygon.

INPUT:

- ``po`` -- a polygon, which can be either of type HPolygon or the refined type HPolygonOpt

OUTPUT:

List of vertices as an array of vertex pairs, Array{Array{Float64,1},1}.
"""
function vertices_list(po::Union{HPolygon, HPolygonOpt})::Array{Array{Float64,1},1}
    n = length(po.constraints)
    vlist = [intersection(Line(po.constraints[i]), Line(po.constraints[i+1])) for i = 1:n-1]
    push!(vlist, intersection(Line(po.constraints[n]), Line(po.constraints[1])))
    return vlist
end


"""
    plot_polygon(P, backend, [name], [gridlines])

Plot a polygon given in constraint form.

INPUT:

- ``P`` -- a polygon, given as a HPolygon or the refined class HPolygonOpt

- ``backend`` -- (optional, default: ``'pyplot'``): select the plot backend; valid
  options are:

                 -  ``'pyplot_savefig'`` -- use PyPlot package, save to a file

                 -  ``'pyplot_inline'`` -- use PyPlot package, showing in external program

                 - ``'gadfly'`` -- use Gadfly package, showing in browser

                 - ``''`` -- (empty string), return nothing, without plotting

- ``name`` -- (optional, default: ``'plot.png'``) the filename of the plot
  (if it is saved to disk)

- ``gridlines`` -- (optional, default: false) to display or not gridlines in
                   the output plot

EXAMPLES:

This function can receive one polygon, as in:

    julia> using LazySets, PyPlot
    julia> H = HPolygon([LinearConstraint([1.0, 0.0], 0.6), LinearConstraint([0.0, 1.0], 0.6),
           LinearConstraint([-1.0, 0.0], -0.4), LinearConstraint([0.0, -1.0], -0.4)])
    julia> plot_polygon(H, backend="pyplot_inline");

Multiple polygons can be plotted passing a list instead of a single element:

    julia> Haux = HPolygon([LinearConstraint([1.0, 0.0], 1.2), LinearConstraint([0.0, 1.0], 1.2),
           LinearConstraint([-1.0, 0.0], -0.8), LinearConstraint([0.0, -1.0], -0.8)])
    julia> plot_polygon([H, Haux], backend="pyplot_inline");
"""
function plot_polygon(P::Union{HPolygon, HPolygonOpt, Array{HPolygon, 1},
        SubArray{HPolygon, 1}, Array{HPolygonOpt, 1}};
        backend="pyplot_savefig", name="plot.png", gridlines=false, color="blue",
        plot_labels::Vector{String}=["", ""])
    if backend in ["pyplot_inline", "pyplot_savefig"]
        if !isdefined(:PyPlot)
            error("this backend requires that your script loads the PyPlot module")
        end
        eval(Expr(:using, :PyPlot))
        if backend == "pyplot_savefig"
            PyPlot.ioff()  # turn off interactive plotting
            fig = PyPlot.figure()
        end
        gridlines ? grid("on") : grid("off")
        # check if P is iterable
        applicable(start, P) ? nothing : P = [P]
        PyPlot.xlabel(plot_labels[1])
        PyPlot.ylabel(plot_labels[2])
        for pi in P
            vlist = hcat(vertices_list(pi)...).'  # each ROW represents a vertex
            # we repeat the 1st element (to "close" the rectangle)
            xcoords = [vlist[:, 1]; vlist[1, 1]]
            ycoords = [vlist[:, 2]; vlist[1, 2]]
            PyPlot.plot(xcoords, ycoords, color=color, linewidth=0.6)
        end
        if backend == "pyplot_savefig"
            PyPlot.savefig(name, bbox_inches="tight")
            PyPlot.close(fig)
        end

    elseif backend in ["gadfly"]
        if !isdefined(:Gadfly)
            error("this backend requires that your script loads the Gadfly module")
        end
        eval(Expr(:using, :Gadfly))
        # check if P is iterable
        applicable(start, P) ? nothing : P = [P]
        layers_array = Vector{Vector{Gadfly.Layer}}(length(P))
        for i in eachindex(P)
            vlist = hcat(vertices_list(P[i])...).'  # each ROW represents a vertex
            layers_array[i] = Gadfly.layer(x=vlist[:, 1], y=vlist[:, 2], Geom.polygon(preserve_order=false, fill=true))
        end
        Gadfly.plot(layers_array...)

    elseif backend == ""
        return nothing
    else
        error("plot backend not valid")
    end

end

export HPolygon, HPolygonOpt, addconstraint!, is_contained, plot_polygon,
       VPolygon, tovrep, vertices_list
