#= concrete implementations of binary intersections between sets =#

export intersection

for T in [:LazySet, :AbstractSingleton, :Interval, :Universe, :LinearMap]
    @eval begin
        function intersection(∅::EmptySet, X::$T)
            @assert dim(∅) == dim(X) "cannot take the intersection between a " *
                "$(dim(∅))-dimensional empty set and a $(dim(X))-dimensional set"
            return ∅
        end

        # symmetric method
        function intersection(X::$T, ∅::EmptySet)
            @assert dim(∅) == dim(X) "cannot take the intersection between a " *
                "$(dim(∅))-dimensional empty set and a $(dim(X))-dimensional set"
            return ∅
        end
    end
end

# disambiguation
function intersection(∅₁::EmptySet, ∅₂::EmptySet)
    @assert dim(∅₁) == dim(∅₂) "cannot take the intersection between two " *
            "empty sets of dimensions $(dim(∅₁)) and $(dim(∅₂))"
    return ∅₁
end

"""
    intersection(S::AbstractSingleton, X::LazySet)

Return the intersection of a singleton with another set.

### Input

- `S` -- singleton
- `X` -- another set

### Output

If the sets intersect, the result is `S`.
Otherwise, the result is the empty set.
"""
@commutative function intersection(S::AbstractSingleton, X::LazySet)
    return _intersection_singleton(S, X)
end

function _intersection_singleton(S::AbstractSingleton, X)
    N = promote_type(eltype(S), eltype(X))
    return element(S) ∈ X ? S : EmptySet{N}(dim(S))
end

# disambiguation
function intersection(S1::AbstractSingleton, S2::AbstractSingleton)
    N = promote_type(eltype(S1), eltype(S1))
    return _isapprox(element(S1), element(S2)) ? S1 : EmptySet{N}(dim(S1))
end

"""
    intersection(L1::Line2D, L2::Line2D)

Return the intersection of two two-dimensional lines.

### Input

- `L1` -- first line
- `L2` -- second line

### Output

Three outcomes are possible:

- If the lines are identical, the result is the first line.
- If the lines are parallel and not identical, the result is the empty set.
- Otherwise the result is the only intersection point.

### Algorithm

We first check whether the lines are parallel.
If not, we use [Cramer's rule](https://en.wikipedia.org/wiki/Cramer%27s_rule)
to compute the intersection point.

### Examples

The line ``y = -x + 1`` intersected with the line ``y = x``:

```jldoctest
julia> intersection(Line2D([-1., 1.], 0.), Line2D([1., 1.], 1.))
Singleton{Float64, Vector{Float64}}([0.5, 0.5])

julia> intersection(Line2D([1., 1.], 1.), Line2D([1., 1.], 1.))
Line2D{Float64, Vector{Float64}}([1.0, 1.0], 1.0)
```
"""
function intersection(L1::Line2D, L2::Line2D)
    det = _det(L1, L2)
    if isapproxzero(det)
        if isapprox(L1.b, L2.b) # lines are identical
            return L1
        else
            N = promote_type(eltype(L1), eltype(L2))
            return EmptySet{N}(dim(L1)) # lines are disjoint
        end
    else # intersection is a point
        @inbounds begin
            det⁻¹ = 1/det
            x = (L1.b * L2.a[2] - L1.a[2] * L2.b) * det⁻¹
            y = (L1.a[1] * L2.b - L1.b * L2.a[1]) * det⁻¹
        end
        return Singleton([x, y])
    end
end

"""
    intersection(a::LineSegment, b::Line2D)

Compute the intersection of a line and a line segment in two dimensions.

### Input

- `a` -- LineSegment
- `b` -- Line2D

### Output

If the sets do not intersect, the result is the empty set.
Otherwise the result is the singleton or line segment that describes the intersection.
"""
function intersection(a::LineSegment, b::Line2D)
    # cast a as line
    ap = Line2D(a.p, a.q)
    # find intersection between a' and b
    m = intersection(ap, b)
    if m == ap
        # if this equals a, then all of the segment is the intersection
        return a
    elseif m isa Singleton && m.element ∈ a
        # if the intersection between lines is in the segment
        return m
    else
        # no intersection
        N = promote_type(eltype(a), eltype(b))
        return EmptySet{N}(2)
    end
end

# symmetric method
intersection(a::Line2D, b::LineSegment) = intersection(b, a)

"""
    intersection(a::LineSegment, b::LineSegment)

Return the intersection of two two-dimensional line segments.

### Input

- `a` -- first line segment
- `b` -- second line segment

### Output

A singleton, line segment or the empty set depending on the result of the intersection.

### Notes

- If the line segments cross, or are parallel and have one point in common,
  that point is returned.

- If the line segments are parallel and have a line segment in common, that
  segment is returned.

- Otherwise, if there is no intersection, an empty set is returned.
"""
function intersection(a::LineSegment, b::LineSegment)

    # cast each segment as a line
    ap = Line2D(a.p, a.q)
    bp = Line2D(b.p, b.q)

    # find intersection between the lines
    m = intersection(ap, bp)
    N = promote_type(eltype(a), eltype(b))
    if m == ap
        # determine which segment is in both
        p1 = max(min(a.p[1], a.q[1]), min(b.p[1], b.q[1]))
        p2 = max(min(a.p[2], a.q[2]), min(b.p[2], b.q[2]))
        q1 = min(max(a.p[1], a.q[1]), max(b.p[1], b.q[1]))
        q2 = min(max(a.p[2], a.q[2]), max(b.p[2], b.q[2]))
        if _isapprox(p1, q1) && _isapprox(p2, q2)
             # edges have a point in common
             return Singleton([p1, p2])

        elseif _leq(p1, q1) && _leq(p2, q2)
             return LineSegment([p1, p2], [q1, q2])

        else
            # no intersection
            return EmptySet{N}(2)
        end

    elseif m isa Singleton && m.element ∈ a && m.element ∈ b
        # if the intersection between lines is in the segments
        return m

    else
        # no intersection
        return EmptySet{N}(2)
    end
end

"""
    intersection(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle)

Return the intersection of two hyperrectangles.

### Input

- `H1` -- first hyperrectangle
- `H2` -- second hyperrectangle

### Output

If the hyperrectangles do not intersect, the result is the empty set.
Otherwise the result is the hyperrectangle that describes the intersection.

### Algorithm

In each isolated direction `i` we compute the rightmost left border and the
leftmost right border of the hyperrectangles.
If these borders contradict, then the intersection is empty.
Otherwise the result uses these borders in each dimension.
"""
function intersection(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle)
    n = dim(H1)
    N = promote_type(eltype(H1), eltype(H2))
    v_high = Vector{N}(undef, n)
    v_low = Vector{N}(undef, n)
    @inbounds for i in 1:n
        v_high[i] = min(high(H1, i), high(H2, i))
        v_low[i] = max(low(H1, i), low(H2, i))
        if v_high[i] < v_low[i]
            return EmptySet{N}(dim(H1))
        end
    end
    return Hyperrectangle(high=v_high, low=v_low)
end

# disambiguations
intersection(S::AbstractSingleton, H::AbstractHyperrectangle) = _intersection_singleton(S, H)
intersection(H::AbstractHyperrectangle, S::AbstractSingleton) = _intersection_singleton(S, H)

"""
    intersection(x::Interval, y::Interval)

Return the intersection of two intervals.

### Input

- `x` -- first interval
- `y` -- second interval

### Output

If the intervals do not intersect, the result is the empty set.
Otherwise the result is the interval that describes the intersection.
"""
function intersection(x::Interval, y::Interval)
    if min(y) > max(x) || min(x) > max(y)
        N = promote_type(eltype(x), eltype(y))
        return EmptySet{N}(1)
    else
        return Interval(max(min(x), min(y)), min(max(x), max(y)))
    end
end

"""
    intersection(X::Interval, hs::HalfSpace)

Compute the intersection of an interval and a half-space.

### Input

- `X`  -- interval
- `hs` -- half-space

### Output

If the sets do not intersect, the result is the empty set.
If the interval is fully contained in the half-space, the result is the original
interval. Otherwise the result is the interval that describes the intersection.

### Algorithm

We first handle the special case that the normal vector `a` of `hs` is close to
zero. Then we distinguish the cases that `hs` is a lower or an upper bound.
"""
function intersection(X::Interval, hs::HalfSpace)
    @assert dim(hs) == 1 "cannot take the intersection between an interval " *
                         "and a $(dim(hs))-dimensional half-space"

    a = hs.a[1]
    b = hs.b
    N = promote_type(eltype(X), eltype(hs))

    if _isapprox(a, zero(N))
        if _geq(b, zero(N))
            # half-space is universal
            return X
        else
            # half-space is empty
            return EmptySet{N}(1)
        end
    end

    empty, lbound, ubound = _intersection_interval_halfspace(min(X), max(X), a, b, N)

    if empty
        return EmptySet{N}(1)
    else
        return Interval(lbound, ubound)
    end
end

function _intersection_interval_halfspace(lo, hi, a, b, N)
    # idea:
    # ax ≤ b  <=>  x ≤ b/a  (if a > 0)
    # ax ≤ b  <=>  x ≥ b/a  (if a < 0)
    b_over_a = b / a

    empty = false
    if a > zero(N)
        # half-space is an upper bound
        # check whether ax ≤ b for x = lo
        if _leq(lo, b_over_a)
            # new upper bound: min(hi, b_over_a)
            if _leq(hi, b_over_a)
                return (empty, lo, hi)
            else
                return (empty, lo, b_over_a)
            end
        end
    else
        # half-space is a lower bound
        # check whether ax ≤ b for x = hi
        if _geq(hi, b_over_a)
            # new lower bound: max(lo, b_over_a)
            if _leq(b_over_a, lo)
                return (empty, lo, hi)
            else
                return (empty, b_over_a, hi)
            end
        end
    end

    # intersection is empty
    return (true, hi, lo)
end

# symmetric method
intersection(hs::HalfSpace, X::Interval) = intersection(X, hs)

"""
    intersection(X::Interval, hp::Hyperplane)

Compute the intersection of an interval and a hyperplane.

### Input

- `X`  -- interval
- `hp` -- hyperplane

### Output

If the sets do not intersect, the result is the empty set.
Otherwise the result is the singleton that describes the intersection.
"""
function intersection(X::Interval, hp::Hyperplane)
    @assert dim(hp) == 1 "cannot take the intersection between an interval " *
                         "and a $(dim(hp))-dimensional hyperplane"

    # a one-dimensional hyperplane is just a point
    p = hp.b / hp.a[1]
    if _leq(min(X), p) && _leq(p, max(X))
        return Singleton([p])
    else
        N = promote_type(eltype(X), eltype(hp))
        return EmptySet{N}(1)
    end
end

# symmetric method
intersection(hp::Hyperplane, X::Interval) = intersection(X, hp)

"""
    intersection(X::Interval, Y::LazySet)

Compute the intersection of an interval and a convex set.

### Input

- `X` -- interval
- `Y` -- convex set

### Output

If the sets do not intersect, the result is the empty set.
Otherwise the result is the interval that describes the intersection, which may
be of type `Singleton` if the intersection is very small.
"""
function intersection(X::Interval, Y::LazySet)
    return _intersection_interval(X, Y)
end

@commutative function _intersection_interval(X::Interval, Y::LazySet)
    N = promote_type(eltype(X), eltype(Y))
    lower = max(min(X), -ρ(N[-1], Y))
    upper = min(max(X), ρ(N[1], Y))
    if _isapprox(lower, upper)
        return Singleton([lower])
    elseif lower < upper
        return Interval(lower, upper)
    else
        return EmptySet{N}(1)
    end
end

# disambiguations
@commutative intersection(X::Interval, H::AbstractHyperrectangle) = _intersection_interval(X, H)
@commutative intersection(X::Interval, S::AbstractSingleton) = _intersection_singleton(S, X)
@commutative intersection(X::Interval, L::LinearMap) = _intersection_interval(X, L)

# special case of an axis-aligned half-space and a hyperrectangular set
function intersection(B::AbstractHyperrectangle, H::HalfSpace{N, <:SingleEntryVector{N}}) where {N}
    n = dim(H)
    a = H.a
    b = H.b
    i = a.i
    ai = a.v

    v_high_i = high(B, i)
    v_low_i = low(B, i)

    # intersect the half-space with the hyperrectangle's interval side
    (empty, lo, hi) = _intersection_interval_halfspace(v_low_i, v_high_i, ai, b, N)

    if empty
        return EmptySet{N}(n)

    else
        v_low′ = copy(low(B))
        v_low′[i] = lo

        v_high′ = copy(high(B))
        v_high′[i] = hi

        return Hyperrectangle(low=v_low′, high=v_high′)
    end
 end

# symmetric method
intersection(H::HalfSpace{N, <:SingleEntryVector{N}}, B::AbstractHyperrectangle) where {N} = intersection(B, H)

# disambiguations
intersection(H::HalfSpace{N, <:SingleEntryVector{N}}, X::Interval) where {N} = _intersection_interval(X, H)
intersection(X::Interval, H::HalfSpace{N, <:SingleEntryVector{N}}) where {N} = _intersection_interval(X, H)
intersection(S::AbstractSingleton, H::HalfSpace) = _intersection_singleton(S, H)
intersection(H::HalfSpace, S::AbstractSingleton) = _intersection_singleton(S, H)
intersection(S::AbstractSingleton, H::HalfSpace{N, <:SingleEntryVector{N}}) where {N} =
    _intersection_singleton(S, H)
intersection(H::HalfSpace{N, <:SingleEntryVector{N}}, S::AbstractSingleton) where {N} =
     _intersection_singleton(S, H)

"""
    intersection(P1::AbstractHPolygon, P2::AbstractHPolygon, [prune]::Bool=true)

Return the intersection of two polygons in constraint representation.

### Input

- `P1`    -- first polygon
- `P2`    -- second polygon
- `prune` -- (optional, default: `true`) flag for removing redundant constraints

### Output

If the polygons do not intersect, the result is the empty set.
Otherwise the result is the polygon that describes the intersection.

### Algorithm

We just combine the constraints of both polygons.
To obtain a linear-time algorithm, we interleave the constraints.
If there are two constraints with the same normal vector, we choose the tighter
one.

Redundancy of constraints is checked with
[`remove_redundant_constraints!(::AbstractHPolygon)`](@ref).
"""
function intersection(P1::AbstractHPolygon, P2::AbstractHPolygon, prune::Bool=true)

    # all constraints of one polygon are processed; now add the other polygon's
    # constraints
    @inline function add_remaining_constraints!(c, i, c1, i1, duplicates)
        c[i+1:length(c)-duplicates] = c1[i1:length(c1)]
        return true
    end

    # choose the constraint of c1; the directions are different
    @inline function choose_first_diff_dir!(c, i, i1, i2, c1, c2, duplicates)
        c[i] = c1[i1]
        if i1 == length(c1)
            return add_remaining_constraints!(c, i, c2, i2, duplicates)
        end
        return false
    end

    # choose the constraint of c1; the directions are equivalent (i.e., linearly
    # dependent)
    @inline function choose_first_same_dir!(c, i, i1, i2, c1, c2, duplicates)
        c[i] = c1[i1]
        if i1 == length(c1)
            if i2 < length(c2)
                return add_remaining_constraints!(c, i, c2, i2+1, duplicates)
            end
            return true
        elseif i2 == length(c2)
            return add_remaining_constraints!(c, i, c1, i1+1, duplicates)
        end
        return false
    end

    c1 = constraints_list(P1)
    c2 = constraints_list(P2)
    if length(c1) == 0
        return P2
    elseif length(c2) == 0
        return P1
    end

    # TODO: use common vector type of P1 and P2, see #2046
    N = promote_type(eltype(P1), eltype(P2))
    c = Vector{LinearConstraint{N, Vector{N}}}(undef, length(c1) + length(c2))
    i1 = 1
    i2 = 1
    duplicates = 0
    for i in 1:length(c)
        if c1[i1].a <= c2[i2].a
            if c2[i2].a <= c1[i1].a
                duplicates += 1
                # constraints have the same normal vector: take the tighter one
                if is_tighter_same_dir_2D(c1[i1], c2[i2])
                    # first constraint is tighter
                    if choose_first_same_dir!(c, i, i1, i2, c1, c2, duplicates)
                        break
                    end
                else
                    # second constraint is tighter
                    if choose_first_same_dir!(c, i, i2, i1, c2, c1, duplicates)
                        break
                    end
                end
                i1 += 1
                i2 += 1
            else
                # first constraint comes first
                if choose_first_diff_dir!(c, i, i1, i2, c1, c2, duplicates)
                    break
                end
                i1 += 1
            end
        else
            # second constraint comes first
            if choose_first_diff_dir!(c, i, i2, i1, c2, c1, duplicates)
                break
            end
            i2 += 1
        end
    end
    if duplicates > 0
        deleteat!(c, length(c)-duplicates+1:length(c))
    end

    # TODO see #2187: the above code does *not* sort the constraints correctly.
    # after fixing the above code, we should pass sort_constraints=false again
    P = HPolygon(c, sort_constraints=true)
    if prune
        remove_redundant_constraints!(P)
        if isempty(P)
            return EmptySet{N}(2)
        end
    end
    return P
end

using MathProgBase.SolverInterface: AbstractMathProgSolver

"""
    intersection(P1::AbstractPolyhedron{N},
                 P2::AbstractPolyhedron{N};
                 [backend]=default_lp_solver(N)) where {N}

Compute the intersection of two polyhedra.

### Input

- `P1`      -- polyhedron
- `P2`      -- polyhedron
- `backend` -- (optional, default: `default_lp_solver(N)`) the LP solver used
               for the removal of redundant constraints; see the `Notes` section
               below for details

### Output

An `HPolyhedron` resulting from the intersection of `P1` and `P2`, with the
redundant constraints removed, or an empty set if the intersection is empty.
If one of the arguments is a polytope, the result is an `HPolytope` instead.

### Notes

The default value of the solver backend is `default_lp_solver(N)` and it is used
to run a feasiblity LP to remove the redundant constraints of the intersection.

If you want to use the `Polyhedra` library, pass an appropriate backend. For
example, to use the default Polyhedra library use
`default_polyhedra_backend(P)` or use `CDDLib.Library()` for the CDD library.

There are some shortcomings of the removal of constraints using the default
Polyhedra library; see e.g. #1038 and Polyhedra#146. It is safer to check for
emptiness of intersection before calling this function in those cases.

### Algorithm

This implementation unifies the constraints of the two sets obtained from the
`constraints_list` method.
"""
intersection(P1::AbstractPolyhedron{N},
             P2::AbstractPolyhedron{N};
             backend=default_lp_solver(N)) where {N} = _intersection_poly(P1, P2, backend=backend)

function _intersection_poly(P1::AbstractPolyhedron{N},
                            P2::AbstractPolyhedron{N};
                            backend=default_lp_solver(N)) where {N}

    # if one of P1 or P2 is bounded => the result is bounded
    HPOLY = (P1 isa AbstractPolytope || P2 isa AbstractPolytope) ?
        HPolytope : HPolyhedron

    # concatenate the linear constraints
    clist_P1 = _normal_Vector(P1) # TODO fix to similar type
    clist_P2 = _normal_Vector(P2)
    Q = HPOLY([clist_P1; clist_P2])

    # remove redundant constraints
    if backend isa AbstractMathProgSolver
        # if Q is empty => the feasiblity LP for the list of constraints of Q
        # is infeasible and remove_redundant_constraints! returns `false`
        if remove_redundant_constraints!(Q, backend=backend)
            return Q
        else
            return EmptySet{N}(dim(P1))
        end
    else
        # the correct way for this condition would be to check if `backend`
        # isa Polyhedra.Library; since that would require using Polyhedra: Library
        # and it is an optional dependency we opt to fallback without checking

        # convert to a Polyhedra's hrep
        Qph = polyhedron(Q; backend=backend)

        # remove the redundancies
        removehredundancy!(Qph)

        if isempty(Qph)
            return EmptySet{N}(dim(P1))
        else
            # convert back to HPOLY
            return convert(HPOLY, Qph)
        end
    end
end

# disambiguations
intersection(S::AbstractSingleton, P::AbstractPolyhedron) = _intersection_singleton(S, P)
intersection(P::AbstractPolyhedron, S::AbstractSingleton) = _intersection_singleton(S, P)
intersection(X::Interval, P::AbstractPolyhedron) = _intersection_interval(X, P)
intersection(P::AbstractPolyhedron, X::Interval) = _intersection_interval(X, P)

"""
    intersection(P1::Union{VPolygon, VPolytope}, P2::Union{VPolygon, VPolytope};
                 [backend]=nothing,
                 [prunefunc]=removevredundancy!)

Compute the intersection of two polytopes in vertex representation.

### Input

- `P1`        -- polytope in vertex representation
- `P2`        -- polytope in vertex representation
- `backend`   -- (optional, default: `nothing`) the backend for polyhedral
                 computations
- `prunefunc` -- (optional, default: `removevredundancy!`) function to prune
                 the vertices of the result

### Output

A `VPolytope`.
"""
function intersection(P1::Union{VPolygon, VPolytope},
                      P2::Union{VPolygon, VPolytope};
                      backend=nothing,
                      prunefunc=nothing)
    n = dim(P1)
    @assert n == dim(P2) "expected polytopes with equal dimensions but they " *
                         "are $(dim(P1)) and $(dim(P2)) respectively"

    # fast path for one and two-dimensional sets
    if n == 1
        Q1 = overapproximate(P1, Interval)
        Q2 = overapproximate(P2, Interval)
        Pint = intersection(Q1, Q2)
        return convert(VPolytope, Pint)
    elseif n == 2
        v1 = convex_hull(vertices_list(P1))
        v2 = convex_hull(vertices_list(P2))
        v12 = _intersection_vrep_2d(v1, v2)
        return VPolytope(v12)
    end

    if isnothing(backend)
        backend = default_polyhedra_backend(P1)
    end
    if isnothing(prunefunc)
        prunefunc = removevredundancy!
    end

    # general case: convert to half-space representation
    Q1 = polyhedron(P1; backend=backend)
    Q2 = polyhedron(P2; backend=backend)
    Pint = Polyhedra.intersect(Q1, Q2)
    prunefunc(Pint)
    return VPolytope(Pint)
end

"""
    intersection(P1::VPolygon, P2::VPolygon; apply_convex_hull::Bool=true)

Compute the intersection of two polygons in vertex representation.

### Input

- `P1` -- polygon in vertex representation
- `P2` -- polygon in vertex representation
- `apply_convex_hull` -- (default, optional: `true`) use the flag to skip the
                         computation of the convex hull in the resulting `VPolygon`

### Output

A `VPolygon` or an `EmptySet` if the intersection is empty.

### Algorithm

This function applies the [Sutherland–Hodgman polygon
clipping algorithm](https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm).
The implementation is based on the one found in
[rosetta code](http://www.rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Julia).
"""
function intersection(P1::VPolygon, P2::VPolygon; apply_convex_hull::Bool=true)
    v1 = vertices_list(P1)
    v2 = vertices_list(P2)
    v12 = _intersection_vrep_2d(v1, v2)

    if isempty(v12)
        N = promote_type(eltype(P1), eltype(P2))
        return EmptySet{N}(2)
    else
        return VPolygon(v12, apply_convex_hull=apply_convex_hull)
    end
end

"""
    intersection(cup::UnionSet, X::LazySet)

Return the intersection of a union of two convex sets and another convex set.

### Input

- `cup` -- union of two convex sets
- `X`   -- convex set

### Output

The union of the pairwise intersections, expressed as a `UnionSet`.
If one of those sets is empty, only the other set is returned.
"""
function intersection(cup::UnionSet, X::LazySet)
    return intersection(cup.X, X) ∪ intersection(cup.Y, X)
end

# symmetric method
function intersection(X::LazySet, cup::UnionSet)
    return intersection(cup, X)
end

# disambiguation
intersection(S::AbstractSingleton, cup::UnionSet) = _intersection_singleton(S, cup)
intersection(cup::UnionSet, S::AbstractSingleton) = _intersection_singleton(S, cup)

"""
    intersection(cup::UnionSetArray, X::LazySet)

Return the intersection of a union of a finite number of convex sets and another
convex set.

### Input

- `cup` -- union of a finite number of convex sets
- `X`   -- convex set

### Output

The union of the pairwise intersections, expressed as a `UnionSetArray`.
"""
function intersection(cup::UnionSetArray, X::LazySet)
    return UnionSetArray([intersection(Y, X) for Y in array(cup)])
end

# symmetric method
intersection(X::LazySet, cup::UnionSetArray) = intersection(cup, X)

# disambiguation
intersection(S::AbstractSingleton, cup::UnionSetArray) = _intersection_singleton(S, cup)
intersection(cup::UnionSetArray, S::AbstractSingleton) = _intersection_singleton(S, cup)

"""
    intersection(L::LinearMap, S::LazySet)

Return the intersection of a lazy linear map and a convex set.

### Input

- `L` -- linear map
- `S` -- convex set

### Output

The polytope obtained by the intersection of `l.M * L.X` and `S`.
"""
function intersection(L::LinearMap, S::LazySet)
    return intersection(linear_map(L.M, L.X), S)
end

# symmetric method
intersection(S::LazySet, L::LinearMap) = intersection(L, S)

# disambiguation
function intersection(L1::LinearMap, L2::LinearMap)
    return intersection(linear_map(L1.M, L1.X), linear_map(L2.M, L2.X))
end
intersection(S::AbstractSingleton, L::LinearMap) = _intersection_singleton(S, L)
intersection(L::LinearMap, S::AbstractSingleton) = _intersection_singleton(S, L)

"""
    intersection(U::Universe, X::LazySet)

Return the intersection of a universe and a convex set.

### Input

- `U` -- universe
- `X` -- convex set

### Output

The set `X`.
"""
function intersection(U::Universe, X::LazySet)
    return X
end

# symmetric method
intersection(X::LazySet, U::Universe) = X

# disambiguations
intersection(U::Universe, ::Universe) = U
intersection(U::Universe, P::AbstractPolyhedron) = P
intersection(P::AbstractPolyhedron, U::Universe) = P
intersection(U::Universe, S::AbstractSingleton) = S
intersection(S::AbstractSingleton, U::Universe) = S
intersection(X::Interval, U::Universe) = X
intersection(U::Universe, X::Interval) = X
intersection(L::LinearMap, U::Universe) = L
intersection(U::Universe, L::LinearMap) = L

"""
    intersection(P::AbstractPolyhedron, rm::ResetMap)

Return the intersection of a polyhedron and a polyhedral reset map.

### Input

- `P`  -- polyhedron
- `rm` -- polyhedral reset map

### Output

A polyhedron.

### Notes

We assume that `rm` is polyhedral, i.e., has a `constraints_list` method
defined.
"""
function intersection(P::AbstractPolyhedron, rm::ResetMap)
    return intersection(P, HPolyhedron(constraints_list(rm)))
end

# symmetric method
intersection(rm::ResetMap, P::AbstractPolyhedron) = intersection(P, rm)

# more efficient version for polytopic
function intersection(P::AbstractPolyhedron,
                      rm::ResetMap{N, <:AbstractPolytope}) where {N}
    return intersection(P, HPolytope(constraints_list(rm)))
end

# symmetric method
function intersection(rm::ResetMap{N, <:AbstractPolytope},
                      P::AbstractPolyhedron) where {N}
    return intersection(P, rm)
end

intersection(U::Universe, X::CartesianProductArray) = X

# symmetric method
intersection(X::CartesianProductArray, U::Universe) = intersection(U, X)

# disambiguation
intersection(P::AbstractSingleton, rm::ResetMap) = _intersection_singleton(S, rm)
intersection(rm::ResetMap, P::AbstractSingleton) = _intersection_singleton(S, rm)
intersection(P::AbstractSingleton, rm::ResetMap{N, <:AbstractPolytope}) where {N} =
    _intersection_singleton(S, rm)
intersection(rm::ResetMap{N, <:AbstractPolytope}, P::AbstractSingleton) where {N} =
    _intersection_singleton(S, rm)
intersection(U::Universe, rm::ResetMap) = rm
intersection(rm::ResetMap, U::Universe) = rm
intersection(U::Universe, rm::ResetMap{N, <:AbstractPolytope}) where {N} = rm
intersection(rm::ResetMap{N, <:AbstractPolytope}, U::Universe) where {N} = rm
intersection(rm::ResetMap, X::Interval) = _intersection_interval(X, Y)
intersection(X::Interval, rm::ResetMap) = _intersection_interval(X, Y)
intersection(rm::ResetMap{N, <:AbstractPolytope}, X::Interval) where {N} = _intersection_interval(X, Y)
intersection(X::Interval, rm::ResetMap{N, <:AbstractPolytope}) where {N} = _intersection_interval(X, Y)

"""
        intersection(X::CartesianProductArray, Y::CartesianProductArray)

Return the intersection between cartesian products of a finite number of convex sets.

### Input

 - `X` -- cartesian product of a finite number of convex sets
 - `Y` -- cartesian product of a finite number of convex sets

### Output

The decomposed set which represents concrete intersection between `X` and `Y`

### Algorithm

This algorithm intersect corresponding blocks between sets.
"""
function intersection(X::CartesianProductArray, Y::CartesianProductArray)
    @assert same_block_structure(array(X), array(Y)) "block structure has to be the same"

    return CartesianProductArray([intersection(array(X)[i], array(Y)[i]) for i in eachindex(array(X))])
end

"""
    intersection(cpa::CartesianProductArray, P::AbstractPolyhedron)

Compute the intersection of a Cartesian product of a finite number of polyhedral
sets with a polyhedron.

### Input

- `cpa` -- Cartesian product of a finite number of polyhedral sets
- `P`   -- polyhedron

### Output

A Cartesian product of a finite number of polyhedral sets.
See the *Algorithm* section below for details about the structure.

### Notes

The restriction to polyhedral sets in `cpa` only applies to the blocks that are
actually intersected with `P` (see the *Algorithm* section below for details).
All other blocks are not considered by the intersection and remain identical.

### Algorithm

The underlying idea of the algorithm is to exploit the unconstrained dimensions
of `P`.
Without loss of generality, assume that `cpa` has the structure ``X × Y × Z`` such that only the
dimensions of ``Y`` are constrained in ``P``, and denoting a suitable projection
of ``P`` to the dimensions of ``Y`` with ``P|_Y``, we have the following
equivalence:

```math
    (X × Y × Z) ∩ P = X × (Y ∩ P|_Y) × Z
```

Note that ``Y`` may still consist of many blocks.
However, due to the structural restriction of a Cartesian product, we cannot
break down this set further even if ``P|_Y`` is still unconstrained in some
dimensions of blocks in ``Y``.
This would require a restructuring of the dimensions.
Consider this example:

```math
    Y := [0, 1] × [1, 2] × [2, 3]
    P|_Y := x₁ + x₃ ≤ 2
    Y ∩ P|_Y = 0 ≤ x₁ ∧ 1 ≤ x₂ ≤ 2 ∧ 2 ≤ x₃ ∧ x₁ + x₃ ≤ 2
```

Even though the constraints of dimension ``x₂`` are decoupled from the rest,
due to the last constraint the Cartesian product cannot be broken down further.
In particular, the result ``Y ∩ P|_Y`` is a polyhedron in this implementation.

Now we explain the implementation of the above idea.
We first identify the dimensions in which `P` is constrained.
Then we identify the block dimensions of ``X × Y × Z`` such that ``Y`` has
minimal dimension.
Finally, we convert ``Y`` to a polyhedron and intersect it with a suitable
projection of `P`.
"""
function intersection(cpa::CartesianProductArray, P::AbstractPolyhedron)
    return _intersection_cpa_polyhedron(cpa, P)
end

function _intersection_cpa_polyhedron(cpa::CartesianProductArray,
                                      P::AbstractPolyhedron)
    # search for the indices of the block trisection into
    # "unconstrained | constrained | unconstrained" (the first and third section
    # may be empty)
    constrained_dims = constrained_dimensions(P)
    if isempty(constrained_dims)
        return cpa
    end
    minimal = minimum(constrained_dims)
    maximal = maximum(constrained_dims)
    lower_bound = 1
    blocks = array(cpa)
    cb_start = 0
    cb_end = length(blocks)
    dim_start = 0
    dim_end = dim(P)
    for (i, block) in enumerate(blocks)
        dim_block = dim(block)
        upper_bound = lower_bound + dim_block - 1
        interval = lower_bound:upper_bound
        if cb_start == 0 && minimal ∈ interval
            dim_start = lower_bound
            cb_start = i
        end
        if maximal ∈ interval
            dim_end = upper_bound
            cb_end = i
            break
        end
        lower_bound = upper_bound + 1
    end
    # compute intersection with constrained blocks
    cap = _intersection_cpa_polyhedron_constrained(
        CartesianProductArray(blocks[cb_start:cb_end]), P, dim_start:dim_end)

    # construct result
    result_array = vcat(blocks[1:(cb_start - 1)],  # unconstrained
                        [cap],  # constrained, intersected
                        blocks[(cb_end + 1):length(blocks)])  # unconstrained
    return CartesianProductArray(result_array)
end

function _intersection_cpa_polyhedron_constrained(cpa::CartesianProductArray,
                                                  P::AbstractPolyhedron,
                                                  constrained_dims)
    T = isbounded(cpa) ? HPolytope : HPolyhedron
    hpoly_low_dim = T(constraints_list(cpa))
    cap_low_dim = intersection(hpoly_low_dim, project(P, constrained_dims))
end

# symmetric method
intersection(P::AbstractPolyhedron, cpa::CartesianProductArray) = _intersection_cpa_polyhedron(cpa, P)

# disambiguation
intersection(S::AbstractSingleton, cpa::CartesianProductArray) = _intersection_singleton(S, cpa)
intersection(cpa::CartesianProductArray, S::AbstractSingleton) = _intersection_singleton(S, cpa)
intersection(cpa::CartesianProductArray, X::Interval) = _intersection_interval(X, cpa)
intersection(X::Interval, cpa::CartesianProductArray) = _intersection_interval(X, cpa)

"""
        intersection(Z::AbstractZonotope{N}, H::HalfSpace{N}; backend=default_lp_solver(N)) where {N}

Return the intersection between a zonotopic set and a halfspace.

### Input

 - `Z` -- zonotopic set
 - `H` -- halfspace

### Output

If the sets do not intersect, the output is the empty set, if the zonotopic set
is fully contained in the halfspace, the zonotopic set is returned, otherwise the
output is the concrete intersection between `Z` and `H`.

### Algorithm

First there is a disjointness test, if that is false, there is an inclusion test,
if that is false then the concrete intersection is computed.
"""
function intersection(Z::AbstractZonotope{N}, H::HalfSpace{N}; backend=default_lp_solver(N)) where {N}
    n = dim(Z)
    isdisjoint(Z, H) && return EmptySet(n)
    issubset(Z, H) && return Z
    return _intersection_poly(Z, H, backend=backend)
end

# symmetric method
intersection(H::HalfSpace{N}, Z::AbstractZonotope{N}; backend=default_lp_solver(N)) where {N} = intersection(Z, H, backend=backend)
