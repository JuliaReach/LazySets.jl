"""
    overapproximate(S::LazySet{N}, ::Type{<:HPolygon})::HPolygon where {N<:Real}

Return an approximation of a given 2D set as a box-shaped polygon.

### Input

- `S` -- set, assumed to be two-dimensional
- `HPolygon` for dispatch

### Output

A box-shaped polygon in constraint representation.
"""
function overapproximate(S::LazySet{N},
                         ::Type{<:HPolygon})::HPolygon where {N<:Real}
    @assert dim(S) == 2
    pe, pn, pw, ps = box_bounds(S)
    constraints = Vector{LinearConstraint{eltype(pe)}}(4)
    constraints[1] = LinearConstraint(DIR_EAST(N), dot(pe, DIR_EAST(N)))
    constraints[2] = LinearConstraint(DIR_NORTH(N), dot(pn, DIR_NORTH(N)))
    constraints[3] = LinearConstraint(DIR_WEST(N), dot(pw, DIR_WEST(N)))
    constraints[4] = LinearConstraint(DIR_SOUTH(N), dot(ps, DIR_SOUTH(N)))
    return HPolygon(constraints)
end

"""
    overapproximate(S::LazySet)::HPolygon

Alias for `overapproximate(S, HPolygon)`.
"""
overapproximate(S::LazySet)::HPolygon = overapproximate(S, HPolygon)

"""
    overapproximate(S::LazySet, Type{<:Hyperrectangle})::Hyperrectangle

Return an approximation of a given 2D set as a hyperrectangle.

### Input

- `S` -- set, assumed to be two-dimensional
- `Hyperrectangle` for dispatch

### Output

A hyperrectangle.
"""
function overapproximate(S::LazySet, ::Type{<:Hyperrectangle})::Hyperrectangle
    @assert dim(S) == 2
    pe, pn, pw, ps = box_bounds(S)
    radius = [(pe[1] - pw[1]) / 2, (pn[2] - ps[2]) / 2]
    center = [pw[1] + radius[1], ps[2] + radius[2]]
    return Hyperrectangle(center, radius)
end

# helper function
@inline function box_bounds(S::LazySet{N}) where {N<:Real}
    # evaluate support vector on box directions
    pe = σ(DIR_EAST(N), S)
    pn = σ(DIR_NORTH(N), S)
    pw = σ(DIR_WEST(N), S)
    ps = σ(DIR_SOUTH(N), S)
    return (pe, pn, pw, ps)
end

"""
    overapproximate(S::LazySet, ɛ::Real)::HPolygon

Return an ɛ-close approximation of the given 2D set (in terms of Hausdorff
distance) as a polygon.

### Input

- `S` -- set, assumed to be two-dimensional
- `ɛ` -- error bound

### Output

A polygon in constraint representation.
"""
function overapproximate(S::LazySet, ɛ::Real)::HPolygon
    @assert dim(S) == 2

    return tohrep(approximate(S, ɛ))
end

"""
    overapproximate(S::ConvexHull{N, Zonotope{N}, Zonotope{N}},
                     ::Type{<:Zonotope})::Zonotope where {N<:Real}

Overapproximate the convex hull of two zonotopes.

### Input

- `S` -- convex set, assumed to be two-dimensional
- `Zonotope` for dispatch

### Algorithm

This function implements the method proposed in
*Reachability of Uncertain Linear Systems Using Zonotopes, A. Girard, HSCC 2005*.

It should be noted that the output zonotope is not necessarily the minimal enclosing
zonotope, which is in general expensive in high dimensions. This is further investigated
in: *Zonotopes as bounding volumes, L. J. Guibas et al, Proc. of Symposium on Discrete Algorithms, pp. 803-812*.
"""
function overapproximate(S::ConvexHull{N, Zonotope{N}, Zonotope{N}},
                         ::Type{<:Zonotope})::Zonotope where {N<:Real}
    Z1, Z2 = S.X, S.Y
    center = (Z1.center+Z2.center)/2
    generators = [(Z1.generators .+ Z2.generators) (Z1.center - Z2.center) (Z1.generators .- Z2.generators)]/2
    return Zonotope(center, generators)
end
