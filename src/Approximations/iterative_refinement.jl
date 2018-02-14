"""
    LocalApproximation{N<:Real}
Type that represents a local approximation in 2D.
### Fields
- `p1`         -- first inner point
- `d1`         -- first direction
- `p2`         -- second inner point
- `d2`         -- second direction
- `q`          -- intersection of the lines l1 ⟂ d1 at p1 and l2 ⟂ d2 at p2
- `refinable`  -- states if this approximation is refinable
- `err`        -- error upper bound
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
- `approx_list` -- vector of local approximations
"""
struct PolygonalOverapproximation{N<:Real}
    S::LazySet
    approx_list::Vector{LocalApproximation{N}}
end

# initialize a polygonal overapproximation with an empty list
PolygonalOverapproximation{N}(S::LazySet) where {N<:Real} =
    PolygonalOverapproximation(S::LazySet, LocalApproximation{N}[])

"""
    new_approx(S::LazySet, p1::Vector{N}, d1::Vector{N}, p2::Vector{N}, d2::Vector{N}) where {N<:Real}
### Input
- `S`          -- convex set
- `p1`         -- first inner point
- `d1`         -- first direction
- `p2`         -- second inner point
- `d2`         -- second direction
### Output
A local approximation of `S` in the given directions.
"""
function new_approx(S::LazySet, p1::Vector{N}, d1::Vector{N}, p2::Vector{N}, d2::Vector{N}) where {N<:Real}
    if norm(p1-p2, 2) < TOL(N)
        # this approximation cannot be refined and we set q = p1 by convention
        ap = LocalApproximation{N}(p1, d1, p2, d2, p1, false, zero(N))
    else
        ndir = normalize([p2[2]-p1[2], p1[1]-p2[1]])
        q = intersection(Line(d1, dot(d1, p1)), Line(d2, dot(d2, p2)))
        approx_error = min(norm(q - σ(ndir, S)), dot(ndir, q - p1))
        refinable = (approx_error > TOL(N)) && !(norm(p1-q, 2) < TOL(N) || norm(q-p2, 2) < TOL(N))
        ap = LocalApproximation{N}(p1, d1, p2, d2, q, refinable, approx_error)
    end
    return ap
end

"""
    addapproximation!(Ω::PolygonalOverapproximation,
        p1::Vector{N}, d1::Vector{N}, p2::Vector{N}, d2::Vector{N}) where {N<:Real}
### Input
- `Ω`          -- polygonal overapproximation of a convex set
- `p1`         -- first inner point
- `d1`         -- first direction
- `p2`         -- second inner point
- `d2`         -- second direction
### Output
The list of local approximations in `Ω` of the set `Ω.S` is updated in-place and
the new approximation is returned by this function.
"""
function addapproximation!(Ω::PolygonalOverapproximation,
    p1::Vector{N}, d1::Vector{N}, p2::Vector{N}, d2::Vector{N})::LocalApproximation{N} where {N<:Real}

    approx = new_approx(Ω.S, p1, d1, p2, d2)
    push!(Ω.approx_list, approx)
    return approx
end

"""
    refine(Ω::PolygonalOverapproximation, i::Int)::Tuple{LocalApproximation, LocalApproximation}
Refine a given local approximation of the polygonal approximation of a convex set,
by splitting along the normal direction to the approximation.
### Input
- `Ω`   -- polygonal overapproximation of a convex set
- `i`   -- integer index for the local approximation to be refined
### Output
The tuple consisting of the refined right and left local approximations.
"""
function refine(Ω::PolygonalOverapproximation, i::Int)::Tuple{LocalApproximation, LocalApproximation}
    R = Ω.approx_list[i]
    @assert R.refinable

    ndir = normalize([R.p2[2]-R.p1[2], R.p1[1]-R.p2[1]])
    s = σ(ndir, Ω.S)
    ap1 = new_approx(Ω.S, R.p1, R.d1, s, ndir)
    ap2 = new_approx(Ω.S, s, ndir, R.p2, R.d2)
    return (ap1, ap2)
end

"""
    tohrep(Ω::PolygonalOverapproximation{N})::AbstractHPolygon where {N<:Real}
Convert a polygonal overapproximation into a concrete polygon.
### Input
- `Ω`   -- polygonal overapproximation of a convex set
### Output
A polygon in constraint representation.
"""
function tohrep(Ω::PolygonalOverapproximation{N})::AbstractHPolygon where {N<:Real}
    p = HPolygon{N}()
    for ai in Ω.approx_list
        addconstraint!(p, LinearConstraint(ai.d1, dot(ai.d1, ai.p1)))
    end
    return p
end

"""
    approximate(S::LazySet{N},
                ɛ::N)::PolygonalOverapproximation{N} where {N<:Real}
Return an ɛ-close approximation of the given 2D convex set (in terms of
Hausdorff distance) as an inner and an outer approximation composed by sorted
local `Approximation2D`.
### Input
- `S` -- 2D convex set
- `ɛ` -- error bound
### Output
An ɛ-close approximation of the given 2D convex set.
"""
function approximate(S::LazySet{N},
                     ɛ::N)::PolygonalOverapproximation{N} where {N<:Real}

    # initialize box directions
    pe = σ(DIR_EAST(N), S)
    pn = σ(DIR_NORTH(N), S)
    pw = σ(DIR_WEST(N), S)
    ps = σ(DIR_SOUTH(N), S)

    Ω = PolygonalOverapproximation{N}(S)

    addapproximation!(Ω, pe, DIR_EAST(N), pn, DIR_NORTH(N))
    addapproximation!(Ω, pn, DIR_NORTH(N), pw, DIR_WEST(N))
    addapproximation!(Ω, pw, DIR_WEST(N), ps, DIR_SOUTH(N))
    addapproximation!(Ω, ps, DIR_SOUTH(N), pe, DIR_EAST(N))

    i = 1
    while i <= length(Ω.approx_list)
        if Ω.approx_list[i].err <= ɛ
            # if this approximation doesn't need to be refined, consider the next
            # one in the queue (counter-clockwise order wrt d1)
            # if the approximation is not refinable => continue
            i += 1
        else
            inext = i + 1

            (la1, la2) = refine(Ω, i)

            Ω.approx_list[i] = la1

            redundant = inext > length(Ω.approx_list) ? false : (norm(la2.p1-Ω.approx_list[inext].p1) < TOL(N)) && (norm(la2.q-Ω.approx_list[inext].q) < TOL(N))
            if redundant
                # if it is redundant, keep the refined approximation
                Ω.approx_list[inext] = la2
            else
                insert!(Ω.approx_list, inext, la2)
            end

        end
    end
    return Ω
end
