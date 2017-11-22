"""
Type that represents a local approximation in 2D.

### Fields

- `p1`        -- the first inner point
- `d1`        -- the first direction
- `p2`        -- the second inner point
- `d2`        -- the second direction
- `err`       -- the error made
- `ndir`      -- a normal direction of the inner approximation
- `refinable` -- states if this approximation is refinable
"""
struct Approximation2D{N<:AbstractFloat}
    p1::Vector{N}
    d1::Vector{N}
    p2::Vector{N}
    d2::Vector{N}
    err::N
    ndir::Vector{N}
    refinable::Bool
end

function Approximation2D(p1::Vector{N}, d1::Vector{N}, p2::Vector{N}, d2::Vector{N}) where{N<:AbstractFloat}
    ndir = [p2[2]-p1[2], p1[1]-p2[1]]
    norm_ndir = norm(ndir)

    if norm_ndir > TOL_DIR
        ndir = ndir/norm_ndir
        q = intersection(Line(d1, dot(d1, p1)), Line(d2, dot(d2, p2)))
        Approximation2D(p1, d1, p2, d2, dot(ndir, q) - dot(ndir, p1), ndir, true)
    else
        Approximation2D(p1, d1, p2, d2, zero(N), ndir, false)
    end
end

"""
    refine(X, A)

Refine the given approximation.

### Input

- `X`      -- set which is approximated
- `approx` -- approximation to refine
"""
function refine(X::LazySet, approx::Approximation2D)
    q = σ(approx.ndir, X)
    (Approximation2D(approx.p1, approx.d1, q, approx.ndir), Approximation2D(q, approx.ndir, approx.p2, approx.d2))
end

"""
    approximate(X, ɛ)

Return an ɛ-close approximation of the given 2D set (in terms of Hausdorff
distance) as an inner and an outer approximation composed by sorted local
`Approximation2D`.

### Input

- `X` -- a 2D set defined by its support function
- `ɛ` -- the error bound
"""
function approximate(X::LazySet, ɛ::Float64)::Vector{Approximation2D}

    # start with box directions
    pe = σ(DIR_EAST, X)
    pn = σ(DIR_NORTH, X)
    pw = σ(DIR_WEST, X)
    ps = σ(DIR_SOUTH, X)
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
            (la1, la2) = refine(X, queue[i])
            queue[i] = la1
            insert!(queue, i+1, la2)
        end
    end
    return queue
end
