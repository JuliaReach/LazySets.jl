export intersection!

@commutative function intersection(∅::EmptySet, X::LazySet)
    return _intersection_emptyset(∅, X)
end

function intersection(S1::AbstractSingleton, S2::AbstractSingleton)
    N = promote_type(eltype(S1), eltype(S1))
    return _isapprox(element(S1), element(S2)) ? S1 : EmptySet{N}(dim(S1))
end

@commutative function intersection(S::AbstractSingleton, H::AbstractHyperrectangle)
    return _intersection_singleton(S, H)
end

@commutative function intersection(S::AbstractSingleton, X::LazySet)
    return _intersection_singleton(S, X)
end

function _intersection_singleton(S::AbstractSingleton, X)
    @assert dim(S) == dim(X) "cannot take the intersection between a " *
                             "$(dim(S))-dimensional set and a $(dim(X))-dimensional set"
    N = promote_type(eltype(S), eltype(X))
    return element(S) ∈ X ? S : EmptySet{N}(dim(S))
end

# this method can also be called with `HalfSpace` arguments
function _intersection_line2d(L1, L2)
    det = right_turn(L1.a, L2.a)
    if isapproxzero(det)
        if isapprox(L1.b, L2.b) # lines are identical
            return _to_Line2D(L1)
        else
            N = promote_type(eltype(L1), eltype(L2))
            return EmptySet{N}(dim(L1)) # lines are disjoint
        end
    else # intersection is a point
        @inbounds begin
            x = (L1.b * L2.a[2] - L1.a[2] * L2.b) / det
            y = (L1.a[1] * L2.b - L1.b * L2.a[1]) / det
        end
        return Singleton([x, y])
    end
end

_to_Line2D(L::Line2D) = L
_to_Line2D(H::HalfSpace) = Line2D(H.a, H.b)

@commutative function intersection(LS::LineSegment, L2::Line2D)
    # cast LS as line
    L1 = Line2D(LS.p, LS.q)
    # find intersection between the lines
    m = intersection(L1, L2)
    if m == L1
        # if the lines are equal, then the intersection is the whole line segment
        return LS
    elseif m isa Singleton && m.element ∈ LS
        # the intersection between the lines is in the segment
        return m
    else
        # no intersection
        N = promote_type(eltype(LS), eltype(L2))
        return EmptySet{N}(2)
    end
end

"""
    intersection(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle)

Compute the intersection of two hyperrectangular sets.

### Input

- `H1` -- hyperrectangular set
- `H2` -- hyperrectangular set

### Output

If the hyperrectangular sets do not intersect, the result is the empty set.
Otherwise the result is the hyperrectangle that describes the intersection.

### Algorithm

In each isolated direction `i` we compute the rightmost left border and the
leftmost right border of the hyperrectangular sets.
If these borders contradict, then the intersection is empty.
Otherwise the resulting hyperrectangle uses these borders in each dimension.
"""
function intersection(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle)
    n = dim(H1)
    @assert n == dim(H2) "cannot take the intersection between a " *
                         "$n-dimensional set and a $(dim(H2))-dimensional set"
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
    return Hyperrectangle(; high=v_high, low=v_low)
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
@commutative function intersection(X::Interval, hs::HalfSpace)
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
    # idea: assuming that a != 0, we have
    # ax ≤ b  <=>  x ≤ b/a  (if a > 0)
    # ax ≤ b  <=>  x ≥ b/a  (if a < 0)
    b_over_a = b / a

    empty = false
    if a > zero(N)
        # half-space is an upper bound
        # check whether ax ≤ b for x = lo
        if _leq(lo, b_over_a)
            # new upper bound: min(hi, b_over_a)
            if hi < b_over_a
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
            if b_over_a < lo
                return (empty, lo, hi)
            else
                return (empty, b_over_a, hi)
            end
        end
    end

    # intersection is empty
    return (true, hi, lo)
end

@commutative function intersection(X::Interval, hp::Hyperplane)
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
@commutative function intersection(X::Interval, Y::LazySet)
    return _intersection_interval(X, Y)
end

function _intersection_interval(X::Interval, Y::LazySet)
    @assert dim(Y) == 1 "cannot take the intersection between an interval " *
                        "and a $(dim(Y))-dimensional set"
    @assert isconvextype(typeof(Y)) "this implementation requires a convex " *
                                    "set, but got $(typeof(Y))"

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

# special case of an axis-aligned half-space and a hyperrectangular set
@commutative function intersection(B::AbstractHyperrectangle,
                                   H::HalfSpace{<:Any,<:SingleEntryVector})
    n = dim(H)
    a = H.a
    b = H.b
    i = a.i
    ai = a.v
    N = promote_type(eltype(B), eltype(H))

    if _isapprox(ai, zero(N))
        if _geq(b, zero(N))
            # half-space is universal
            return X
        else
            # half-space is empty
            return EmptySet{N}(dim(H))
        end
    end

    v_high_i = high(B, i)
    v_low_i = low(B, i)

    # intersect the half-space with the hyperrectangle's interval side
    empty, lo, hi = _intersection_interval_halfspace(v_low_i, v_high_i, ai, b, N)

    if empty
        return EmptySet{N}(n)

    else
        v_low′ = copy(low(B))
        v_low′[i] = lo

        v_high′ = copy(high(B))
        v_high′[i] = hi

        return Hyperrectangle(; low=v_low′, high=v_high′)
    end
end

"""
    intersection(P1::AbstractHPolygon, P2::AbstractHPolygon; [prune]::Bool=true)

Compute the intersection of two polygons in constraint representation.

### Input

- `P1`    -- polygon in constraint representation
- `P2`    -- polygon in constraint representation
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
function intersection(P1::AbstractHPolygon, P2::AbstractHPolygon;
                      prune::Bool=true)

    # all constraints of one polygon are processed; now add the other polygon's
    # constraints
    @inline function add_remaining_constraints!(c, i, c1, i1, duplicates)
        c[(i + 1):(length(c) - duplicates)] = c1[i1:length(c1)]
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
                return add_remaining_constraints!(c, i, c2, i2 + 1, duplicates)
            end
            return true
        elseif i2 == length(c2)
            return add_remaining_constraints!(c, i, c1, i1 + 1, duplicates)
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
    c = Vector{HalfSpace{N,Vector{N}}}(undef, length(c1) + length(c2))
    i1 = 1
    i2 = 1
    duplicates = 0
    for i in eachindex(c)
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
        deleteat!(c, (length(c) - duplicates + 1):length(c))
    end

    # TODO see #2187: the above code does *not* sort the constraints correctly.
    # after fixing the above code, we should pass sort_constraints=false again
    P = HPolygon(c; sort_constraints=true)
    if prune
        remove_redundant_constraints!(P)
        if isempty(P)
            return EmptySet{N}(2)
        end
    end
    return P
end

# use H-rep algorithm for mixed polygons
function intersection(P::AbstractPolygon, Q::AbstractPolygon; prune::Bool=true)
    return intersection(tohrep(P), tohrep(Q); prune=prune)
end

"""
    intersection(P1::AbstractPolyhedron{N}, P2::AbstractPolyhedron{N};
                 [backend]=default_lp_solver(N), [prune]::Bool=true) where {N}

Compute the intersection of two polyhedra.

### Input

- `P1`      -- polyhedron
- `P2`      -- polyhedron
- `backend` -- (optional, default: `default_lp_solver(N)`) the LP solver used
               for the removal of redundant constraints; see the *Notes* section
               below for details
- `prune`   -- (optional, default: `true`) flag for removing redundant
               constraints

### Output

An `HPolyhedron` resulting from the intersection of `P1` and `P2`, with the
redundant constraints removed, or an empty set if the intersection is empty.
If one of the arguments is a polytope, the result is an `HPolytope` instead.

### Notes

The default value of the solver backend is `default_lp_solver(N)` and it is used
to run a feasiblity LP to remove the redundant constraints of the intersection.

If you want to use the `Polyhedra` library, pass an appropriate backend. For
example, use `default_polyhedra_backend(P)` for the default Polyhedra library,
or use `CDDLib.Library()` for the CDD library.

There are some shortcomings of the removal of constraints using the default
Polyhedra library; see e.g. #1038 and Polyhedra#146. It is safer to check for
emptiness of intersection before calling this function in those cases.

### Algorithm

This implementation unifies the constraints of the two sets obtained from the
`constraints_list` method.
"""
function intersection(P1::AbstractPolyhedron{N},
                      P2::AbstractPolyhedron{N};
                      backend=default_lp_solver(N),
                      prune::Bool=true) where {N}
    return _intersection_poly(P1, P2; backend=backend, prune=prune)
end

function _intersection_poly(P1::AbstractPolyhedron{N},
                            P2::AbstractPolyhedron{N};
                            backend=default_lp_solver(N),
                            prune::Bool=true) where {N}

    # if one of P1 or P2 is bounded => the result is bounded
    HPOLY = (P1 isa AbstractPolytope || P2 isa AbstractPolytope) ?
            HPolytope : HPolyhedron

    # concatenate the linear constraints
    clist_P1 = _normal_Vector(P1) # TODO fix to similar type
    clist_P2 = _normal_Vector(P2)
    Q = HPOLY([clist_P1; clist_P2])

    # remove redundant constraints
    if _is_polyhedra_backend(backend)
        # convert to Polyhedra's hrep
        Qph = polyhedron(Q; backend=backend)

        # remove the redundancies
        if prune
            removehredundancy!(Qph)
        end

        if isempty(Qph)
            return EmptySet{N}(dim(P1))
        else
            # convert back to HPOLY
            return convert(HPOLY, Qph)
        end
    else
        # if Q is empty => the feasiblity LP for the list of constraints of Q
        # is infeasible and remove_redundant_constraints! returns `false`
        if !prune || remove_redundant_constraints!(Q; backend=backend)
            return Q
        else
            return EmptySet{N}(dim(P1))
        end
    end
end

"""
    intersection(P1::Union{VPolygon, VPolytope}, P2::Union{VPolygon, VPolytope};
                 [backend]=nothing, [prunefunc]=nothing)

Compute the intersection of two polytopes in vertex representation.

### Input

- `P1`        -- polytope in vertex representation
- `P2`        -- polytope in vertex representation
- `backend`   -- (optional, default: `nothing`) the backend for polyhedral
                 computations
- `prunefunc` -- (optional, default: `nothing`) function to prune the vertices
                 of the result

### Output

A `VPolytope`.

### Notes

If `prunefunc` is `nothing`, this implementation sets it to
`(X -> removevredundancy!(X; tol=_ztol(eltype(P1))))`.
"""
function intersection(P1::Union{VPolygon,VPolytope},
                      P2::Union{VPolygon,VPolytope};
                      backend=nothing, prunefunc=nothing)
    n = dim(P1)
    @assert n == dim(P2) "expected polytopes with equal dimensions but they " *
                         "are $(dim(P1)) and $(dim(P2)) respectively"

    # fast path for one- and two-dimensional sets
    if n == 1
        Q1 = overapproximate(P1, Interval)
        Q2 = overapproximate(P2, Interval)
        Pint = intersection(Q1, Q2)
        return isempty(Pint) ? Pint : convert(VPolytope, Pint)
    elseif n == 2
        v1 = convex_hull(vertices_list(P1))
        v2 = convex_hull(vertices_list(P2))
        v12 = _intersection_vrep_2d(v1, v2)
        return isempty(v12) ? v12 : VPolytope(v12)
    end

    if isnothing(backend)
        backend = default_polyhedra_backend(P1)
    end

    # general case: convert to half-space representation
    Q1 = polyhedron(P1; backend=backend)
    Q2 = polyhedron(P2; backend=backend)
    Pint = Polyhedra.intersect(Q1, Q2)

    N = promote_type(eltype(P1), eltype(P2))
    if isnothing(prunefunc)
        prunefunc = (X -> _removevredundancy!(X; N=N))
    end
    prunefunc(Pint)

    if isempty(Pint)
        return EmptySet{N}(n)
    end
    return convert(VPolytope, Pint)
end

intersection(cup1::UnionSet, cup2::UnionSet) = _intersection_us(cup1, cup2)

"""
# Extended help

    intersection(cup::UnionSet, X::LazySet)

### Output

The union of the pairwise intersections, expressed as a `UnionSet`.
If one of those sets is empty, only the other set is returned.
"""
@commutative function intersection(cup::UnionSet, X::LazySet)
    return _intersection_us(cup, X)
end

function _intersection_us(cup::UnionSet, X::LazySet)
    return intersection(first(cup), X) ∪ intersection(second(cup), X)
end

intersection(cup1::UnionSetArray, cup2::UnionSetArray) = _intersection_usa(cup1, cup2)

@commutative function intersection(cup::UnionSetArray, X::LazySet)
    return _intersection_usa(cup, X)
end

function _intersection_usa(cup::UnionSetArray, X::LazySet)
    sets = [intersection(Y, X) for Y in cup]
    l = length(sets)
    filter!(!isempty, sets)
    if length(sets) > 1
        if length(sets) < l
            sets = [X for X in sets]  # re-allocate set-specific array
        end
        return UnionSetArray(sets)
    elseif length(sets) == 1
        return sets[1]
    else
        N = promote_type(eltype(cup), eltype(X))
        return EmptySet{N}(dim(X))
    end
end

function intersection(L1::LinearMap, L2::LinearMap)
    return intersection(linear_map(L1.M, L1.X), linear_map(L2.M, L2.X))
end

@commutative function intersection(L::LinearMap, X::LazySet)
    return intersection(linear_map(L.M, L.X), X)
end

@commutative function intersection(U::Universe, X::LazySet)
    return _intersection_universe(U, X)
end

"""
# Extended help

    intersection(P::AbstractPolyhedron, rm::ResetMap)

### Notes

This method assumes that `rm` is polyhedral, i.e., has a `constraints_list`
method defined.
"""
@commutative function intersection(P::AbstractPolyhedron, rm::ResetMap)
    return intersection(P, HPolyhedron(constraints_list(rm)))
end

# more specific version for polytopic reset map
@commutative function intersection(P::AbstractPolyhedron,
                                   rm::ResetMap{N,<:AbstractPolytope}) where {N}
    return intersection(P, HPolytope(constraints_list(rm)))
end

"""
        intersection(X::CartesianProductArray, Y::CartesianProductArray)

Compute the intersection between Cartesian products of a finite number of sets
with identical decomposition.

### Input

 - `X` -- Cartesian product of a finite number of sets
 - `Y` -- Cartesian product of a finite number of sets

### Output

The decomposed set that represents the concrete intersection of `X` and `Y`.

### Algorithm

This algorithm intersects the corresponding blocks of the sets.
"""
function intersection(X::CartesianProductArray, Y::CartesianProductArray)
    @assert same_block_structure(array(X), array(Y)) "the block structure " *
                                                     "has to be identical"

    return CartesianProductArray([intersection(array(X)[i], array(Y)[i])
                                  for i in eachindex(array(X))])
end

"""
    intersection(cpa::Union{CartesianProduct,CartesianProductArray},
                 P::AbstractPolyhedron)

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
of `P`. Without loss of generality, assume that `cpa` has the structure
``X × Y × Z`` such that only the dimensions of ``Y`` are constrained in ``P``.
By denoting a suitable projection of ``P`` to the dimensions of ``Y`` with
``P|_Y``, we have the following equivalence:

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
due to the last constraint, the Cartesian product cannot be broken down further.
In particular, the result ``Y ∩ P|_Y`` is a polyhedron in this implementation.

Now we explain the implementation of the above idea.
We first identify the dimensions in which `P` is constrained.
Then we identify the block dimensions of ``X × Y × Z`` such that ``Y`` has
minimal dimension.
Finally, we convert ``Y`` to a polyhedron and intersect it with a suitable
projection of `P`.
"""
@commutative function intersection(cpa::CartesianProductArray, P::AbstractPolyhedron)
    return _intersection_cpa(cpa, P)
end

@commutative function intersection(cpa::CartesianProduct, P::AbstractPolyhedron)
    return _intersection_cpa(cpa, P)
end

function _intersection_cpa(cpa::Union{CartesianProduct,CartesianProductArray},
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
    cb_start = 0
    cb_end = length(cpa)
    dim_start = 0
    dim_end = dim(P)
    for (i, block) in enumerate(cpa)
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
    cap = _intersection_polyhedron_constrained(CartesianProductArray(cpa[cb_start:cb_end]), P,
                                               dim_start:dim_end)

    # construct result
    result_array = vcat(cpa[1:(cb_start - 1)],  # unconstrained
                        [cap],  # constrained, intersected
                        cpa[(cb_end + 1):end])  # unconstrained
    return CartesianProductArray(result_array)
end

function _intersection_polyhedron_constrained(X::LazySet, P::AbstractPolyhedron,
                                              constrained_dims)
    T = isboundedtype(typeof(X)) ? HPolytope : HPolyhedron
    hpoly_low_dim = T(constraints_list(X))
    return cap_low_dim = intersection(hpoly_low_dim, project(P, constrained_dims))
end

"""
    intersection(Z::AbstractZonotope{N}, H::HalfSpace{N};
                 [backend]=default_lp_solver(N), [prune]::Bool=true) where {N}

Compute the intersection between a zonotopic set and a half-space.

### Input

- `Z`       -- zonotopic set
- `H`       -- half-space
- `backend` -- (optional, default: `default_lp_solver(N)`) the LP solver used
               for the removal of redundant constraints
- `prune`   -- (optional, default: `true`) flag for removing redundant
               constraints

### Output

If the sets do not intersect, the output is the empty set. If the zonotopic set
is fully contained in the half-space, the zonotopic set is returned. Otherwise,
the output is the concrete intersection between `Z` and `H`.

### Algorithm

First there is a disjointness test. If that is negative, there is an inclusion
test. If that is negative, then the concrete intersection is computed.
"""
@commutative function intersection(Z::AbstractZonotope{N}, H::HalfSpace{N};
                                   backend=default_lp_solver(N),
                                   prune::Bool=true) where {N}
    n = dim(Z)
    isdisjoint(Z, H) && return EmptySet{N}(n)
    issubset(Z, H) && return Z
    return _intersection_poly(Z, H; backend=backend, prune=prune)
end

"""
    intersection!(X::Star, H::HalfSpace)

Compute the intersection between a star set and a half-space, in-place.

### Input

- `X` -- star set
- `H` -- half-space

### Output

The modified star set.
"""
function intersection!(X::Star, H::HalfSpace)
    _intersection_star!(center(X), basis(X), predicate(X), H)
    return X
end

function _intersection_star!(c, V, P::Union{HPoly,HPolygon,HPolygonOpt}, H::HalfSpace)
    a′ = transpose(V) * H.a
    b′ = H.b - dot(H.a, c)
    H′ = HalfSpace(a′, b′)
    return addconstraint!(P, H′)
end

@commutative function intersection(X::Star, H::HalfSpace)
    c = center(X)
    V = basis(X)
    N = eltype(X)
    Pnew = convert(HPolyhedron{N,Vector{N}}, predicate(X))
    Xnew = Star(c, V, Pnew)
    return intersection!(Xnew, H)
end

# sets in H-rep copy the set instead
@commutative function intersection(X::Star{N,VN,MN,PT},
                                   H::HalfSpace) where {N,
                                                        VN<:AbstractVector{N},MN<:AbstractMatrix{N},
                                                        PT<:Union{HPoly,HPolygon,HPolygonOpt}}
    return intersection!(copy(X), H)
end

"""
    _bound_intersect_2D(Z::Zonotope, L::Line2D)

Evaluate the support function in the direction [0, 1] of the intersection
between the given zonotope and line.

### Input

- `Z` -- zonotope
- `L` -- vertical 2D line

### Output

The support function in the direction [0, 1] of the intersection between the
given zonotope and line.

### Notes

The algorithm assumes that the given line is vertical and that the intersection
between the given sets is not empty.

### Algorithm

This function implements [LeGuernic09; Algorithm 8.2](@citet).
"""
function _bound_intersect_2D(Z::Zonotope, L::Line2D)
    c = center(Z)
    P = copy(c)
    G = genmat(Z)
    r = ngens(Z)
    g(x) = view(G, :, x)
    for i in 1:r
        gi = g(i)
        if !isupwards(gi)
            gi .= -gi
        end
        P .= P - gi
    end
    G = sortslices(G; dims=2, by=x -> atan(x[2], x[1])) # sort gens
    if P[1] < L.b
        G .= G[:, end:-1:1]
    end
    j = 1
    while isdisjoint(LineSegment(P, P + 2g(j)), L)
        P .= P + 2g(j)
        j += 1
        if j > size(G, 2)
            error("got an unexpected error; check that the sets intersect")
        end
    end
    singleton = intersection(LineSegment(P, P + 2g(j)), L)
    return element(singleton)[2]
end

# ============== #
# disambiguation #
# ============== #

for T in (:AbstractSingleton, :Interval, :Universe, :LinearMap, :UnionSet, :UnionSetArray)
    @eval @commutative function intersection(∅::EmptySet, X::$T)
        return _intersection_emptyset(∅, X)
    end
end

for T in (:AbstractPolyhedron, :AbstractSingleton, :Interval, :LinearMap, :ResetMap,
          :(ResetMap{<:Any,<:AbstractPolytope}), :CartesianProduct, :CartesianProductArray,
          :UnionSet, :UnionSetArray)
    @eval @commutative function intersection(U::Universe, X::$T)
        return _intersection_universe(U, X)
    end
end

for T in (:AbstractPolyhedron, :AbstractHyperrectangle, :(HalfSpace{<:Any,<:SingleEntryVector}),
          :LinearMap, :ResetMap, :(ResetMap{<:Any,<:AbstractPolytope}), :CartesianProduct,
          :CartesianProductArray, :UnionSet)
    @eval @commutative function intersection(X::Interval, Y::($T))
        return _intersection_interval(X, Y)
    end
end

for T in (:AbstractPolyhedron, :Interval, :HalfSpace, :(HalfSpace{<:Any,<:SingleEntryVector}),
          :LinearMap, :ResetMap, :(ResetMap{<:Any,<:AbstractPolytope}), :CartesianProduct,
          :CartesianProductArray, :UnionSet, :UnionSetArray)
    @eval @commutative function intersection(S::AbstractSingleton, X::($T))
        return _intersection_singleton(S, X)
    end
end

for T in (:LinearMap, :UnionSetArray)
    @eval @commutative function intersection(cup::UnionSet, X::($T))
        return _intersection_us(cup, X)
    end
end

for T in (:Interval, :LinearMap)
    @eval @commutative function intersection(cup::UnionSetArray, X::($T))
        return _intersection_usa(cup, X)
    end
end
