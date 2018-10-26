"""
    LocalApproximation{N<:Real}

Type that represents a local approximation in 2D.

### Fields

- `p1`        -- first inner point
- `d1`        -- first direction
- `p2`        -- second inner point
- `d2`        -- second direction
- `q`         -- intersection of the lines l1 ⟂ d1 at p1 and l2 ⟂ d2 at p2
- `refinable` -- states if this approximation is refinable
- `err`       -- error upper bound

### Notes

The criteria for refinable are determined in the method `new_approx`.
"""
struct LocalApproximation{N<:Real}
    p1::Vector{N}
    d1::Vector{N}
    p2::Vector{N}
    d2::Vector{N}
    q::Vector{N}
    refinable::Bool
    err::N
end

"""
    PolygonalOverapproximation{N<:Real}

Type that represents the polygonal approximation of a convex set.

### Fields

- `S`           -- convex set
- `approx_list` -- vector of local approximations that are not finished yet
- `constraints` -- vector of linear constraints that are already finished
"""
struct PolygonalOverapproximation{N<:Real}
    S::LazySet{N}
    approx_list::Vector{LocalApproximation{N}}
    constraints::Vector{LinearConstraint{N}}

    PolygonalOverapproximation(S::LazySet{N}) where N<:Real = new{N}(
        S, Vector{LocalApproximation{N}}(), Vector{LinearConstraint{N}}())
end

"""
    new_approx(S::LazySet, p1::Vector{N}, d1::Vector{N}, p2::Vector{N},
               d2::Vector{N}) where N<:Real

### Input

- `S`  -- convex set
- `p1` -- first inner point
- `d1` -- first direction
- `p2` -- second inner point
- `d2` -- second direction

### Output

A local approximation of `S` in the given directions.
"""
function new_approx(S::LazySet, p1::Vector{N}, d1::Vector{N}, p2::Vector{N},
                    d2::Vector{N}) where N<:Real
    if norm(p1-p2, 2) <= TOL(N)
        # this approximation cannot be refined and we set q = p1 by convention
        ap = LocalApproximation{N}(p1, d1, p2, d2, p1, false, zero(N))
    else
        ndir = normalize([p2[2]-p1[2], p1[1]-p2[1]])
        q = element(intersection(Line(d1, dot(d1, p1)), Line(d2, dot(d2, p2))))
        approx_error = min(norm(q - σ(ndir, S)), dot(ndir, q - p1))
        refinable = (approx_error > TOL(N)) && (norm(p1-q, 2) > TOL(N)) &&
                    (norm(q-p2, 2) > TOL(N))
        ap = LocalApproximation{N}(p1, d1, p2, d2, q, refinable, approx_error)
    end
    return ap
end

"""
    addapproximation!(Ω::PolygonalOverapproximation, p1::Vector{N},
        d1::Vector{N}, p2::Vector{N}, d2::Vector{N}) where N<:Real

### Input

- `Ω`  -- polygonal overapproximation of a convex set
- `p1` -- first inner point
- `d1` -- first direction
- `p2` -- second inner point
- `d2` -- second direction

### Output

The list of local approximations in `Ω` of the set `Ω.S` is updated in-place and
the new approximation is returned by this function.
"""
function addapproximation!(Ω::PolygonalOverapproximation,
                           p1::Vector{N}, d1::Vector{N}, p2::Vector{N},
                           d2::Vector{N})::LocalApproximation{N} where N<:Real

    approx = new_approx(Ω.S, p1, d1, p2, d2)
    push!(Ω.approx_list, approx)
    return approx
end

"""
    refine(Ω::PolygonalOverapproximation,
           R::LocalApproximation
          )::Tuple{LocalApproximation, LocalApproximation}

Refine a given local approximation of the polygonal approximation of a convex
set by splitting along the normal direction to the approximation.

### Input

- `Ω` -- polygonal overapproximation of a convex set
- `i` -- integer index for the local approximation to be refined

### Output

The tuple consisting of the refined right and left local approximations.
"""
function refine(Ω::PolygonalOverapproximation,
                R::LocalApproximation
               )::Tuple{LocalApproximation, LocalApproximation}
    @assert R.refinable

    ndir = normalize([R.p2[2]-R.p1[2], R.p1[1]-R.p2[1]])
    s = σ(ndir, Ω.S)
    ap1 = new_approx(Ω.S, R.p1, R.d1, s, ndir)
    ap2 = new_approx(Ω.S, s, ndir, R.p2, R.d2)
    return (ap1, ap2)
end

"""
    tohrep(Ω::PolygonalOverapproximation{N})::AbstractHPolygon{N} where N<:Real

Convert a polygonal overapproximation into a concrete polygon.

### Input

- `Ω` -- polygonal overapproximation of a convex set

### Output

A polygon in constraint representation.
"""
function tohrep(Ω::PolygonalOverapproximation{N}
               )::AbstractHPolygon{N} where N<:Real
    return HPolygon{N}(Ω.constraints)
end

"""
    approximate(S::LazySet{N},
                ε::N)::PolygonalOverapproximation{N} where N<:Real

Return an ε-close approximation of the given 2D convex set (in terms of
Hausdorff distance) as an inner and an outer approximation composed by sorted
local `Approximation2D`.

### Input

- `S` -- 2D convex set
- `ε` -- error bound

### Output

An ε-close approximation of the given 2D convex set.
"""
function approximate(S::LazySet{N},
                     ε::N)::PolygonalOverapproximation{N} where N<:Real

    # initialize box directions
    pe = σ(DIR_EAST(N), S)
    pn = σ(DIR_NORTH(N), S)
    pw = σ(DIR_WEST(N), S)
    ps = σ(DIR_SOUTH(N), S)

    Ω = PolygonalOverapproximation(S)

    # add constraints in reverse (i.e., clockwise) order
    addapproximation!(Ω, ps, DIR_SOUTH(N), pe, DIR_EAST(N))
    addapproximation!(Ω, pw, DIR_WEST(N), ps, DIR_SOUTH(N))
    addapproximation!(Ω, pn, DIR_NORTH(N), pw, DIR_WEST(N))
    addapproximation!(Ω, pe, DIR_EAST(N), pn, DIR_NORTH(N))

    approx_list = Ω.approx_list
    while !isempty(approx_list)
        approx = pop!(approx_list)

        if !approx.refinable || approx.err <= ε
            # if the approximation is not refinable => continue
            push!(Ω.constraints,
                  LinearConstraint(approx.d1, dot(approx.d1, approx.p1)))
            continue
        end

        # refine
        (la1, la2) = refine(Ω, approx)

        if isempty(approx_list)
            redundant = false
        else
            # check if the next local approximation became redundant
            next = approx_list[end]
            redundant = (norm(la2.p1-next.p1) <= TOL(N)) &&
                (norm(la2.q-next.q) <= TOL(N))
        end
        if redundant
            # replace redundant old constraint
            approx_list[end] = la2
        else
            # add new constraint
            push!(approx_list, la2)
        end
        push!(approx_list, la1)
    end
    return Ω
end
