import Base.<=

export HPolygon, HPolygonOpt, addconstraint!, is_contained,
       VPolygon, tovrep, vertices_list

"""
    HPolygon <: LazySet

Type that represents a convex polygon (in H-representation).

### Fields

- `constraints` --  an array of linear constraints
"""
mutable struct HPolygon <: LazySet
    constraints::Array{LinearConstraint, 1}
end
HPolygon() = HPolygon([])

function dim(P::HPolygon)::Int64
    2
end

"""
    addconstraint!(p, c)

Add a linear constraint to a polygon keeping the constraints sorted by their
normal directions.

### Input

- `p` -- a polygon
- `c` -- the linear constraint to add
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

### Input

- `d` -- direction
- `p` -- polyhedron in H-representation
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

### Input

- `x` -- vector
- `P` -- polygon

### Output

Return rue iff x ∈ P.
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

### Fields

- `P`   -- polygon
- `ind` -- an index in the list of constraints to begin the search to
           evaluate the support functions.

### Notes

This structure is optimized to evaluate the support function/vector with a large
sequence of directions, which are one to one close.
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

### Input

- `P` -- optimized polyhedron in H-representation
"""
function dim(P::HPolygonOpt)::Int64
    return 2
end

"""
    σ(d, P)

Return the support vector of the optimized polygon in a given direction.

### Input

- `d` -- direction
- `P` -- polyhedron in H-representation
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

### Input

- `x` -- vector
- `P` -- polygon

### Output

Return true iff x ∈ P.
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

### Fields

- `vl` -- the list of vertices
"""
mutable struct VPolygon
    vl::Array{Vector{Float64}, 1}
end
VPolygon() = VPolygon([])

"""
    tovrep(s)

Build a vertex representation of the given polygon.

### Input

- `s` -- a polygon in H-representation, HPolygon. The linear constraints are
         assumed sorted by their normal directions.

### Output

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

### Input

- `po` -- a polygon in H-representation. The linear constraints are
          assumed sorted by their normal directions.

### Output

The same polygon in a vertex representation.
"""
function tovrep(po::HPolygonOpt)
    tovrep(po.p)
end

"""
    vertices_list(po)

Return the list of vertices of a convex polygon.

### Input

- `po` -- a polygon, which can be either of type HPolygon or the refined type HPolygonOpt

### Output

List of vertices as an array of vertex pairs, Array{Array{Float64,1},1}.
"""
function vertices_list(po::Union{HPolygon, HPolygonOpt})::Array{Array{Float64,1},1}
    n = length(po.constraints)
    vlist = [intersection(Line(po.constraints[i]), Line(po.constraints[i+1])) for i = 1:n-1]
    push!(vlist, intersection(Line(po.constraints[n]), Line(po.constraints[1])))
    return vlist
end

"""
    plot_Polygon(P::Union{HPolygon, HPolygonOpt}; ...)

Plot a polygon given in constraint form.
 
### Input

- `P` -- a polygon in constraint representation

### Examples

```julia
julia> using LazySets, Plots
julia> H = HPolygon([LinearConstraint([1.0, 0.0], 0.6), LinearConstraint([0.0, 1.0], 0.6),
                     LinearConstraint([-1.0, 0.0], -0.4), LinearConstraint([0.0, -1.0], -0.4)])
julia> plot(H)
```
"""
@recipe function plot_Polygon(P::Union{HPolygon, HPolygonOpt};
                              color="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]
 
     x, y
end

"""
    plot_Polygon(P::Union{Vector{HPolygon}, Vector{HPolygonOpt}}; ...)

Plot an array of polygons given in constraint form.
 
### Input

- `P` -- an array of polygons in constraint representation

### Examples

```julia
julia> using LazySets, Plots
julia> H1 = HPolygon([LinearConstraint([1.0, 0.0], 0.6), LinearConstraint([0.0, 1.0], 0.6),
                      LinearConstraint([-1.0, 0.0], -0.4), LinearConstraint([0.0, -1.0], -0.4)])
julia> H2 = HPolygon([LinearConstraint([2.0, 0.0], 0.6), LinearConstraint([0.0, 2.0], 0.6),
                      LinearConstraint([-2.0, 0.0], -0.4), LinearConstraint([0.0, -2.0], -0.4)])
julia> plot([H1, H2])
```
"""
@recipe function plot_Polygon(P::Union{Vector{HPolygon}, Vector{HPolygonOpt}};
                              seriescolor="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    for Pi in P
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end
