__precompile__(true)

"""
Module `Approximations.jl` -- polygonal approximation of convex sets through
support vectors.
"""
module Approximations

using LazySets

export decompose, overapproximate, approximate,
       radius_approximation, diameter_approximation, box_approximation,
       ballinf_approximation, symmetric_interval_hull

const TOL_DIR = 1e-6

const DIR_EAST = [1., 0.]
const DIR_NORTH = [0., 1.]
const DIR_WEST = [-1., 0.]
const DIR_SOUTH = [0., -1.]

"""
Type that represents a local approximation in 2D.

FIELDS:

- ``p1`` -- the first inner point

- ``d1`` -- the first direction

- ``p2`` -- the second inner point

- ``d2`` -- the second direction

- ``err`` -- the error made

- ``ndir`` -- a normal direction of the inner approximation

- ``refinable`` -- states if this approximation is refinable
"""
struct Approximation2D
    p1::Vector{Float64}
    d1::Vector{Float64}
    p2::Vector{Float64}
    d2::Vector{Float64}
    err::Float64
    ndir::Vector{Float64}
    refinable::Bool
    function Approximation2D(p1::Vector{Float64}, d1::Vector{Float64}, p2::Vector{Float64}, d2::Vector{Float64})
        ndir = [p2[2]-p1[2], p1[1]-p2[1]]
        norm_ndir = norm(ndir)
        #print(norm_ndir, "\n")
        if norm_ndir > TOL_DIR
            ndir = ndir/norm_ndir
            q = intersection(Line(d1, dot(d1, p1)), Line(d2, dot(d2, p2)))
            new(p1, d1, p2, d2, dot(ndir, q) - dot(ndir, p1), ndir, true)
        else
            new(p1, d1, p2, d2, 0., ndir, false)
        end
    end
end

"""
    refine(S, A)

Refine the given approximation.

INPUT:

- ``S`` -- the set which is approximated

- ``approx`` -- the approximation to refine
"""
function refine(S::LazySet, approx::Approximation2D)
    q = σ(approx.ndir, S)
    (Approximation2D(approx.p1, approx.d1, q, approx.ndir), Approximation2D(q, approx.ndir, approx.p2, approx.d2))
end

"""
    approximate(S, ɛ)

Return an ɛ-close approximation of the given 2D set (in terms of Hausdorff
distance) as an inner and an outer approximation composed by sorted local
Approximation2D.

INPUT:

- ``S`` -- a 2D set defined by its support function
- ``ɛ`` -- the error bound
"""
function approximate(S::LazySet, ɛ::Float64)::Array{Approximation2D, 1}
    # box directions
    pe = σ(DIR_EAST, S)
    pn = σ(DIR_NORTH, S)
    pw = σ(DIR_WEST, S)
    ps = σ(DIR_SOUTH, S)
    queue = Approximation2D[]
    push!(queue, Approximation2D(pe, DIR_EAST, pn, DIR_NORTH))
    push!(queue, Approximation2D(pn, DIR_NORTH, pw, DIR_WEST))
    push!(queue, Approximation2D(pw, DIR_WEST, ps, DIR_SOUTH))
    push!(queue, Approximation2D(ps, DIR_SOUTH, pe, DIR_EAST))
    i = 1
    # iterative refinement
    while i <= length(queue)
        if (queue[i].err <= ɛ)
            i += 1
        else
            #print("refining..", "queue error = ", queue[i].err, "\n")
            (la1, la2) = refine(S, queue[i])
            queue[i] = la1
            insert!(queue, i+1, la2)
        end
    end
    return queue
end

"""
    overapproximate(S, ɛ)

Return an ɛ-close approximation of the given 2D set (in terms of Hausdorff
distance) as a polygon.

INPUT:

- ``S`` -- a 2D set defined by its support function
- ``ɛ`` -- the error bound

OUTPUT:

A polygon in constraint representation.
"""
function overapproximate(S::LazySet, ɛ::Float64)::HPolygon
    constraints = [LinearConstr(la.d1, dot(la.d1, la.p1)) for la in approximate(S, ɛ)]
    return HPolygon(constraints)
end

"""
  overapproximate(S)

Return an approximation of the given 2D set as a polygon, using box directions.

INPUT:

- ``S`` -- a 2D set defined by its support function

OUTPUT:

A polygon in constraint representation.
"""
function overapproximate(S::LazySet)::HPolygon
    constraints = Array{LinearConstr, 1}(4)
    # evaluate support vector on box directions
    pe = σ(DIR_EAST, S)
    pn = σ(DIR_NORTH, S)
    pw = σ(DIR_WEST, S)
    ps = σ(DIR_SOUTH, S)
    constraints[1] = LinearConstr(DIR_EAST, dot(pe, DIR_EAST))
    constraints[2] = LinearConstr(DIR_NORTH, dot(pn, DIR_NORTH))
    constraints[3] = LinearConstr(DIR_WEST, dot(pw, DIR_WEST))
    constraints[4] = LinearConstr(DIR_SOUTH, dot(ps, DIR_SOUTH))
    return HPolygon(constraints)
end

"""
    box_approximation(sf)

Overapproximate a set by a box (hyperrectangle). 

INPUT:

``sf`` -- a set

OUTPUT:

``H`` -- a (tight) hyperrectangle

ALGORITHM:

The center of the hyperrectangle is obtained by averaring the support function
the given set in the canonical directions, and the lengths of the sides can be
recovered from the distance among support functions in the same directions.
"""
function box_approximation(sf::LazySet)::Hyperrectangle
    (n, c, r) = box_approximation_helper(sf)
    return Hyperrectangle(c, r)
end


"""
    box_approximation_symmetric(sf)

Overapproximation of a set by a hyperrectangle which contains the origin.

INPUT:

``sf`` -- a set

OUTPUT:

``H`` -- a symmetric interval around the origin which tightly contains the given set

ALGORITHM:

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(sf::LazySet)::Hyperrectangle
    (n, c, r) = box_approximation_helper(sf)
    return Hyperrectangle(zeros(n), abs.(c) + r)
end
# function alias
symmetric_interval_hull = box_approximation_symmetric


"""
    box_approximation_helper(sf)

Common code of box_approximation and box_approximation_symmetric.

INPUT:

``sf`` -- a set

OUTPUT:

``H`` -- a (tight) hyperrectangle

ALGORITHM:

The center of the hyperrectangle is obtained by averaring the support function
the given set in the canonical directions, and the lengths of the sides can be
recovered from the distance among support functions in the same directions.
"""
@inline function box_approximation_helper(sf::LazySet)
    n = dim(sf)
    c = Vector{Float64}(n)
    r = Vector{Float64}(n)
    dplus = zeros(n)
    dminus = zeros(n)
    @inbounds @simd for i in 1:n
        dplus[i] = 1.0
        dminus[i] = -1.0
        htop = ρ(dplus, sf)
        hbottom = -ρ(dminus, sf)
        dplus[i] = 0.0
        dminus[i] = 0.0
        c[i] = (htop+hbottom)/2.
        r[i] = (htop-hbottom)/2.
    end
    return n, c, r
end


"""
    ballinf_approximation(S)

Overapproximation of a set by a ball in the infinity norm.

INPUT:

``sf`` -- a set

OUTPUT:

``H`` -- a ball in the infinity norm which tightly contains the given set

ALGORITHM:

The center and radius of the box are obtained by evaluating the support function
of the given set along the canonical directions.
"""
function ballinf_approximation(sf::LazySet)::BallInf
    n = dim(sf)
    c = Vector{Float64}(n)
    r = 0.
    dplus(i::Int64) = [zeros(i-1); 1.; zeros(n-i)]
    dminus(i::Int64) = [zeros(i-1); -1.; zeros(n-i)]
    @inbounds for i in 1:n
        htop = ρ(dplus(i), sf)
        hbottom = -ρ(dminus(i), sf)
        c[i] = (htop+hbottom)/2.
        rcur = (htop-hbottom)/2.
        if (rcur > r)
            r = rcur
        end
    end
    return BallInf(c, r)
end

"""
    radius_approximation(sf)

Approximate radius of a given set.

This is an approximation in the infinity norm. The radius of a BallInf of center
c and radius r can be approximated by ‖c‖ + r√n, where n is the dimension of the 
vectorspace.

INPUT:

``sf`` -- set
"""
function radius_approximation(sf::LazySet)::Float64
    b = ballinf_approximation(sf)
    return norm(b.center) + b.radius * sqrt(dim(sf))
end

"""
    diameter_approximation(sf)

Approximate diameter of a given set.

The diameter is bounded by 2*radius. Relies on radius_approximation.

INPUT:

- ``sf`` -- set
"""
function diameter_approximation(sf::LazySet)::Float64
    return 2.*radius_approximation(sf)
end


"""
    decompose(X, ɛi)

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace.

INPUT:

- ``X`` -- set represented by support functions
- ``ɛi`` -- array, error bound for each projection

OUTPUT:

A CartesianProductArray corresponding to the cartesian product of 2x2 polygons.

NOTES:

- It assumes blocks of size 2.
- Different error bounds can be passed to different blocks.
- If the type of the given X is CartesianProductArray, we cannot do simply

```
return CartesianProductArray([overapproximate(X.sfarray[i], ɛi[i]) for i in 1:b])
```
because the cartesian product could contain objects in different dimension, a priori.

ALGORITHM:

1. Project the set X into each partition (assuming two-dimensional subspaces, i.e. partitions of the form [2, 2, ..., 2]),
with M*X, where M is the identity matrix in the block coordinates and zero otherwise.
2. Overapproximate the set with a given error bound, ɛi[i], for i = 1,.., b
3. Return the result as an array of support functions.
"""
function decompose(X::LazySet, ɛi::Vector{Float64})::CartesianProductArray
    n = LazySets.dim(X)
    b = div(n, 2)
    result = Array{HPolygon, 1}(b)
    @inbounds for i in 1:b
        M = sparse([1, 2], [2*i-1, 2*i], [1., 1.], 2, n)
        result[i] = overapproximate(M * X, ɛi[i])
    end
    return CartesianProductArray(result)
end

function decompose(X::LazySet, ɛ::Float64)::CartesianProductArray
    n = LazySets.dim(X)
    b = div(n, 2)
    return decompose(X, [ɛ for i in 1:b])
end

function decompose(X::LazySet)::CartesianProductArray
    n = LazySets.dim(X)
    b = div(n, 2)

    DIR_EAST_b = repeat(DIR_EAST, outer=b)
    DIR_NORTH_b = repeat(DIR_NORTH, outer=b)
    DIR_WEST_b = repeat(DIR_WEST, outer=b)
    DIR_SOUTH_b = repeat(DIR_SOUTH, outer=b)

    pe = reshape(σ(DIR_EAST_b, X), 2, :)
    pe = DIR_EAST.' * pe

    pn = reshape(σ(DIR_NORTH_b, X), 2, :)
    pn = DIR_NORTH.' * pn

    pw = reshape(σ(DIR_WEST_b, X), 2, :)
    pw = DIR_WEST.' * pw

    ps = reshape(σ(DIR_SOUTH_b, X), 2, :)
    ps = DIR_SOUTH.' * ps

    result = Array{HPolygon, 1}(b)
    @inbounds for i in 1:b
        result[i] = HPolygon([LinearConstr(DIR_EAST, pe[i]),
                              LinearConstr(DIR_NORTH, pn[i]),
                              LinearConstr(DIR_WEST, pw[i]),
                              LinearConstr(DIR_SOUTH, ps[i])])
    end
    return CartesianProductArray(result)
end

end # module
