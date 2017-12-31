"""
    Approximation2D{N<:AbstractFloat}

Type that represents a local approximation in 2D.

### Fields

- `p1`        -- first inner point
- `d1`        -- first direction
- `p2`        -- second inner point
- `d2`        -- second direction
- `q`         -- intersection of the lines `x : d1.x = p1` and `x : d2.x = p2`
- `err`       -- error made
- `ndir`      -- normal direction of the inner approximation
- `refinable` -- states if this approximation is refinable
"""
struct Approximation2D{N<:AbstractFloat}
    p1::Vector{N}
    d1::Vector{N}
    p2::Vector{N}
    d2::Vector{N}
    q::Vector{N}
    err::N
    ndir::Vector{N}
    refinable::Bool
end

"""
    Approximation2D(p1::Vector{N}, d1::Vector{N}, p2::Vector{N}, d2::Vector{N}) where {N<:AbstractFloat}

Constructor of `Approximation2D` from two inner points and two directions.

### Input

- `p1`        -- first inner point
- `d1`        -- first direction
- `p2`        -- second inner point
- `d2`        -- second direction

### Output

A new `Approximation2D` instance.
"""
function Approximation2D(p1::Vector{N},
                         d1::Vector{N},
                         p2::Vector{N},
                         d2::Vector{N}) where {N<:AbstractFloat}
    ndir = [p2[2]-p1[2], p1[1]-p2[1]]
    norm_ndir = norm(ndir)

    if norm_ndir > TOL_DIR
        ndir = ndir/norm_ndir
        q = intersection(Line(d1, dot(d1, p1)), Line(d2, dot(d2, p2)))
        Approximation2D(p1, d1, p2, d2, q, dot(ndir, q) - dot(ndir, p1), ndir,
                        true)
    else
        Approximation2D(p1, d1, p2, d2, p1, zero(N), ndir, false)
    end
end

"""
    refine(S::LazySet, approx::Approximation2D)::Tuple{Approximation2D, Approximation2D}

Refine the given approximation.

### Input

- `S`      -- 2D convex set that is approximated
- `approx` -- approximation to refine

### Output

The refined approximation.
"""
function refine(S::LazySet,
                approx::Approximation2D)::Tuple{Approximation2D, Approximation2D}
    q = σ(approx.ndir, S)
    return (Approximation2D(approx.p1, approx.d1, q, approx.ndir),
            Approximation2D(q, approx.ndir, approx.p2, approx.d2))
end

"""
    approximate(S::LazySet, ɛ::Float64)::Vector{Approximation2D}

Return an ɛ-close approximation of the given 2D convex set (in terms of
Hausdorff distance) as an inner and an outer approximation composed by sorted
local `Approximation2D`.

### Input

- `S` -- 2D convex set
- `ɛ` -- error bound

### Output

An ɛ-close approximation of the given 2D convex set.
"""
function approximate(S::LazySet, ɛ::Float64)::Vector{Approximation2D}

    # start with box directions
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
    N = length(queue)
    while i <= length(queue)
        if (queue[i].err <= ɛ)
            i += 1

        else
            (la1, la2) = refine(S, queue[i])

            if (maximum(map(abs, la1.q - la2.q)) < TOL_POINTS)
                # update the constraints if a redundance is detected
                la = queue[i]
                queue[i] = Approximation2D(la1.q, la.d1, la.p2, la.d2)
                iM1 = ifelse(i == 1, N, i-1)
                la = queue[iM1]
                queue[iM1] = Approximation2D(la.p1, la.d1, la1.q, la.d2)

            else
                # insert the new constraints if no redundance is detected
                queue[i] = la1
                insert!(queue, i+1, la2)
            end
        end
        N = length(queue)
    end
    return queue
end
