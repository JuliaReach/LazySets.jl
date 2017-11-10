import Base.<=

# Polygons in constraint representation
export HPolygon, HPolygonOpt, addconstraint!, is_contained

"""
    HPolygon <: LazySet

Type that represents a convex polygon in constraint representation, whose edges are
sorted in counter-clockwise fashion with respect to their normal directions.

### Fields

- `constraints` --  an array of linear constraints

### Note

The `HPolygon` constructor *does not perform* sorting of the given list of eddges.
Use `addconstraint!` to iteratively add and sort the edges.
"""
struct HPolygon <: LazySet
    constraints_list::Vector{LinearConstraint}
end
HPolygon() = HPolygon([])

function dim(P::HPolygon)::Int64
    2
end

"""
    addconstraint!(P, constraint)

Add a linear constraint to a polygon in contraing representation keeping the
constraints sorted by their normal directions.

### Input

- `P`          -- a polygon
- `constraint` -- the linear constraint to add
"""
function addconstraint!(P::HPolygon, constraint::LinearConstraint)
    i = length(P.constraints_list)
    while i > 0 && constraint.a <= P.constraints_list[i].a
        i -= 1
    end
    # here P.constraints_list[i] < constraint
    insert!(P.constraints_list, i+1, constraint)
end

"""
    σ(d, P)

Return the support vector of a polygon in a given direction.
Return the zero vector if there are no constraints (i.e., the universal
polytope).

### Input

- `d` -- direction
- `P` -- polyhedron in H-representation
"""
function σ(d::AbstractVector{Float64}, P::HPolygon)::Vector{Float64}
    n = length(P.constraints_list)
    if n == 0
        error("this polygon is empty")
    end
    i = 1
    while i <= n && P.constraints_list[i].a <= d
        i += 1
    end
    if i == 1 || i == n+1
        intersection(Line(P.constraints_list[1]), Line(P.constraints_list[n]))
    else
        intersection(Line(P.constraints_list[i]), Line(P.constraints_list[i-1]))
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
function is_contained(x::AbstractVector{Float64}, P::HPolygon)::Bool
    if (length(x) != 2)
        false
    else
        res = true
        for c in P.constraints_list
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
    constraints_list::Vector{LinearConstraint}
    ind::Int64

    HPolygonOpt(constraints_list) = new(constraints_list, 1)
    HPolygonOpt(constraints_list, ind) = new(constraints_list, ind)
end
HPolygonOpt(H::HPolygon) = HPolygonOpt(H.constraints_list)

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
function σ(d::AbstractVector{Float64}, p::HPolygonOpt)::Vector{Float64}
    n = length(p.constraints_list)
    cl = p.constraints_list
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
        for c in P.constraints_list
            res = res && dot(c.a, x) <= c.b
        end
        res
    end
end

"""
    tovrep(s)

Build a vertex representation of the given polygon.

### Input

- `s` -- a polygon in H-representation, HPolygon

### Output

The same polygon in a vertex representation, VPolygon.

### Note

The linear constraints of the input HPolygon are assumed to be sorted by their
normal directions in counter-clockwise fashion.
"""
function tovrep(s::HPolygon)
    n = length(s.constraints_list)
    p = Vector{Float64}[]
    for i in 1:n-1
        cur = intersection(Line(s.constraints_list[i]), Line(s.constraints_list[i+1]))
        push!(p, cur)
    end
    cur = intersection(Line(s.constraints_list[n]), Line(s.constraints_list[1]))
    push!(p, cur)
    return VPolygon(p)
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
    vertices_list(P)

Return the list of vertices of a convex polygon in constraint representation.

### Input

- `po` -- a polygon, which can be either of type HPolygon or the refined type HPolygonOpt

### Output

List of vertices as an array of vertex pairs, Vector{Vector{Float64}}.
"""
function vertices_list(P::HPolygon)::Vector{Vector{Float64}}
    n = length(P.constraints_list)
    vlist = [intersection(Line(P.constraints_list[i]), Line(P.constraints_list[i+1])) for i = 1:n-1]
    push!(vlist, intersection(Line(P.constraints_list[n]), Line(P.constraints_list[1])))
    return vlist
end

"""
    plot_polygon(P::Union{HPolygon, HPolygonOpt}; ...)

Plot a polygon given in constraint form.

### Input

- `P` -- a polygon in constraint representation

### Examples

```julia
julia> using LazySets, Plots
julia> P = HPolygon([LinearConstraint([1.0, 0.0], 0.6), LinearConstraint([0.0, 1.0], 0.6),
                     LinearConstraint([-1.0, 0.0], -0.4), LinearConstraint([0.0, -1.0], -0.4)])
julia> plot(P)
```
"""
@recipe function plot_polygon(P::Union{HPolygon, HPolygonOpt};
                              color="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]

     x, y
end

"""
    plot_polygons(P::Union{Vector{HPolygon}, Vector{HPolygonOpt}}; ...)

Plot an array of polygons given in constraint form.

### Input

- `P` -- an array of polygons in constraint representation

### Examples

```julia
julia> using LazySets, Plots
julia> P1 = HPolygon([LinearConstraint([1.0, 0.0], 0.6), LinearConstraint([0.0, 1.0], 0.6),
                      LinearConstraint([-1.0, 0.0], -0.4), LinearConstraint([0.0, -1.0], -0.4)])
julia> P2 = HPolygon([LinearConstraint([2.0, 0.0], 0.6), LinearConstraint([0.0, 2.0], 0.6),
                      LinearConstraint([-2.0, 0.0], -0.4), LinearConstraint([0.0, -2.0], -0.4)])
julia> plot([P1, P2])
```
"""
@recipe function plot_polygons(P::Vector{HPolygon};
                              seriescolor="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    for Pi in P
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end
