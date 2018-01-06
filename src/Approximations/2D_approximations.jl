"""
    overapproximate(S::LazySet, ::Type{HPolygon})::HPolygon

Return an approximation of a given 2D convex set as a box-shaped polygon.

### Input

- `S` -- convex set, assumed to be two-dimensional
- `HPolygon` for dispatch

### Output

A box-shaped polygon in constraint representation.
"""
function overapproximate(S::LazySet, ::Type{HPolygon})::HPolygon
    @assert dim(S) == 2
    pe, pn, pw, ps = box_bounds(S)
    constraints = Vector{LinearConstraint{eltype(pe)}}(4)
    constraints[1] = LinearConstraint(DIR_EAST, dot(pe, DIR_EAST))
    constraints[2] = LinearConstraint(DIR_NORTH, dot(pn, DIR_NORTH))
    constraints[3] = LinearConstraint(DIR_WEST, dot(pw, DIR_WEST))
    constraints[4] = LinearConstraint(DIR_SOUTH, dot(ps, DIR_SOUTH))
    return HPolygon(constraints)
end

"""
    overapproximate(S::LazySet)::HPolygon

Alias for `overapproximate(S, HPolygon)`.
"""
overapproximate(S::LazySet)::HPolygon = overapproximate(S, HPolygon)

"""
    overapproximate(S::LazySet, Type{Hyperrectangle})::Hyperrectangle

Return an approximation of a given 2D convex set as a hyperrectangle.

### Input

- `S` -- convex set, assumed to be two-dimensional
- `Hyperrectangle` for dispatch

### Output

A hyperrectangle.
"""
function overapproximate(S::LazySet, ::Type{Hyperrectangle})::Hyperrectangle
    @assert dim(S) == 2
    pe, pn, pw, ps = box_bounds(S)
    radius = [(pe[1] - pw[1]) / 2, (pn[2] - ps[2]) / 2]
    center = [pw[1] + radius[1], ps[2] + radius[2]]
    return Hyperrectangle(center, radius)
end

# helper function
@inline function box_bounds(S)
    # evaluate support vector on box directions
    pe = σ(DIR_EAST, S)
    pn = σ(DIR_NORTH, S)
    pw = σ(DIR_WEST, S)
    ps = σ(DIR_SOUTH, S)
    return (pe, pn, pw, ps)
end

"""
    overapproximate(S::LazySet, ɛ::Float64)::HPolygon

Return an ɛ-close approximation of the given 2D set (in terms of Hausdorff
distance) as a polygon.

### Input

- `S` -- convex set, assumed to be two-dimensional
- `ɛ` -- error bound

### Output

A polygon in constraint representation.
"""
function overapproximate(S::LazySet, ɛ::Float64)::HPolygon
    @assert dim(S) == 2

    return tohrep(approximate(S, ɛ))
end
