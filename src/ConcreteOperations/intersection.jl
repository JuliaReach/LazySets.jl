#= concrete implementations of binary intersections between sets =#

export intersection

"""
    intersection(S::AbstractSingleton{N},
                 X::LazySet{N}
                ) where {N<:Real}

Return the intersection of a singleton with another set.

### Input

- `S` -- singleton
- `X` -- another set

### Output

If the sets intersect, the result is `S`.
Otherwise, the result is the empty set.
"""
function intersection(S::AbstractSingleton{N},
                      X::LazySet{N}
                     ) where {N<:Real}
    return element(S) ∈ X ? S : EmptySet{N}(dim(S))
end

# symmetric method
function intersection(X::LazySet{N},
                      S::AbstractSingleton{N}
                     ) where {N<:Real}
    return intersection(S, X)
end

# disambiguation
function intersection(S1::AbstractSingleton{N},
                      S2::AbstractSingleton{N}
                     ) where {N<:Real}
    return element(S1) == element(S2) ? S1 : EmptySet{N}(dim(S1))
end

"""
    intersection(L1::Line{N}, L2::Line{N}
                ) where {N<:Real}

Return the intersection of two 2D lines.

### Input

- `L1` -- first line
- `L2` -- second line

### Output

If the lines are identical, the result is the first line.
If the lines are parallel and not identical, the result is the empty set.
Otherwise the result is the only intersection point.

### Algorithm

We first check whether the lines are parallel.
If not, we use [Cramer's rule](https://en.wikipedia.org/wiki/Cramer%27s_rule)
to compute the intersection point.

### Examples

The line ``y = -x + 1`` intersected with the line ``y = x``:

```jldoctest
julia> intersection(Line([-1., 1.], 0.), Line([1., 1.], 1.))
Singleton{Float64,Array{Float64,1}}([0.5, 0.5])

julia> intersection(Line([1., 1.], 1.), Line([1., 1.], 1.))
Line{Float64,Array{Float64,1}}([1.0, 1.0], 1.0)
```
"""
function intersection(L1::Line{N}, L2::Line{N}
                     ) where {N<:Real}
    det = L1.a[1]*L2.a[2] - L1.a[2]*L2.a[1]
    # are the lines parallel?
    if isapproxzero(det)
        # are they the same line?
        if isapprox(L1.b, L2.b)
            return L1
        else
            # lines are parallel but not identical
            return EmptySet{N}(dim(L1))
        end
    else
        x = (L1.b*L2.a[2] - L1.a[2]*L2.b)/det
        y = (L1.a[1]*L2.b - L1.b*L2.a[1])/det
        return Singleton([x, y])
    end
end

"""
    intersection(H1::AbstractHyperrectangle{N},
                 H2::AbstractHyperrectangle{N}
                ) where {N<:Real}

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
function intersection(H1::AbstractHyperrectangle{N},
                      H2::AbstractHyperrectangle{N}
                     ) where {N<:Real}
    n = dim(H1)
    c1 = center(H1)
    c2 = center(H2)
    r1 = radius_hyperrectangle(H1)
    r2 = radius_hyperrectangle(H2)
    high = Vector{N}(undef, n)
    low = Vector{N}(undef, n)
    for i in 1:n
        high1 = c1[i] + r1[i]
        low1 = c1[i] - r1[i]
        high2 = c2[i] + r2[i]
        low2 = c2[i] - r2[i]
        high[i] = min(high1, high2)
        low[i] = max(low1, low2)
        if high[i] < low[i]
            return EmptySet{N}(dim(H1))
        end
    end
    return Hyperrectangle(high=high, low=low)
end

# disambiguation
function intersection(S::AbstractSingleton{N},
                      H::AbstractHyperrectangle{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, H)
end
function intersection(H::AbstractHyperrectangle{N},
                      S::AbstractSingleton{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, H)
end

"""
    intersection(x::Interval{N}, y::Interval{N}
                ) where {N<:Real}

Return the intersection of two intervals.

### Input

- `x` -- first interval
- `y` -- second interval

### Output

If the intervals do not intersect, the result is the empty set.
Otherwise the result is the interval that describes the intersection.
"""
function intersection(x::Interval{N}, y::Interval{N}
                     ) where {N<:Real}
    if min(y) > max(x) || min(x) > max(y)
        return EmptySet{N}(1)
    else
        return Interval(max(min(x), min(y)), min(max(x), max(y)))
    end
end

"""
    intersection(X::Interval{N}, hs::HalfSpace{N}
                ) where {N<:Real}

Compute the intersection of an interval and a half-space.

### Input

- `X`  -- interval
- `hs` -- half-space

### Output

If the sets do not intersect, the result is the empty set.
If the interval is fully contained in the half-space, the result is the original
interval.
Otherwise the result is the interval that describes the intersection.

### Algorithm

We first handle the special case that the normal vector `a` of `hs` is close to
zero.
Then we distinguish the cases that `hs` is a lower or an upper bound.
"""
function intersection(X::Interval{N}, hs::HalfSpace{N}
                     ) where {N<:Real}
    @assert dim(hs) == 1 "cannot take the intersection between an interval " *
                         "and a $(dim(hs))-dimensional half-space"

    a = hs.a[1]
    b = hs.b
    if _isapprox(a, zero(N))
        if _geq(b, zero(N))
            # half-space is universal
            return X
        else
            # half-space is empty
            return EmptySet{N}(1)
        end
    end

    # idea:
    # ax ≤ b  <=>  x ≤ b/a  (if a > 0)
    # ax ≤ b  <=>  x ≥ b/a  (if a < 0)
    b_over_a = b / a
    lo = min(X)
    hi = max(X)

    if a > zero(N)
        # half-space is an upper bound
        # check whether ax ≤ b for x = lo
        if _leq(lo, b_over_a)
            # new upper bound: min(hi, b_over_a)
            if _leq(hi, b_over_a)
                return X
            else
                return Interval(lo, b_over_a)
            end
        else
            return EmptySet{N}(1)
        end
    else
        # half-space is a lower bound
        # check whether ax ≤ b for x = hi
        if _geq(hi, b_over_a)
            # new lower bound: max(lo, b_over_a)
            if _leq(b_over_a, lo)
                return X
            else
                return Interval(b_over_a, hi)
            end
        else
            return EmptySet{N}(1)
        end
    end
end

# symmetric method
function intersection(hs::HalfSpace{N}, X::Interval{N}
                     ) where {N<:Real}
    return intersection(X, hs)
end

"""
    intersection(X::Interval{N}, hp::Hyperplane{N}
                ) where {N<:Real}

Compute the intersection of an interval and a hyperplane.

### Input

- `X`  -- interval
- `hp` -- hyperplane

### Output

If the sets do not intersect, the result is the empty set.
Otherwise the result is the singleton that describes the intersection.
"""
function intersection(X::Interval{N}, hp::Hyperplane{N}
                     ) where {N<:Real}
    @assert dim(hp) == 1 "cannot take the intersection between an interval " *
                         "and a $(dim(hp))-dimensional hyperplane"

    # a one-dimensional hyperplane is just a point
    p = hp.b / hp.a[1]
    if _leq(min(X), p) && _leq(p, max(X))
        return Singleton([p])
    else
        return EmptySet{N}(1)
    end
end

# symmetric method
function intersection(hp::Hyperplane{N}, X::Interval{N}
                     ) where {N<:Real}
    return intersection(X, hp)
end

"""
    intersection(X::Interval{N}, Y::LazySet{N}
                ) where {N<:Real}

Compute the intersection of an interval and a convex set.

### Input

- `X` -- interval
- `Y` -- convex set

### Output

If the sets do not intersect, the result is the empty set.
Otherwise the result is the interval that describes the intersection, which may
be of type `Singleton` if the intersection is very small.
"""
function intersection(X::Interval{N}, Y::LazySet{N}
                     ) where {N<:Real}
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

# symmetric method
function intersection(Y::LazySet{N}, X::Interval{N}
                     ) where {N<:Real}
    return intersection(X, Y)
end

# disambiguation
function intersection(X::Interval{N}, H::AbstractHyperrectangle{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{Interval{N}, LazySet{N}}, X, H)
end
function intersection(H::AbstractHyperrectangle{N}, X::Interval{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{Interval{N}, LazySet{N}}, X, H)
end
function intersection(X::Interval{N}, S::AbstractSingleton{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, X)
end
function intersection(S::AbstractSingleton{N}, X::Interval{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, X)
end

"""
    intersection(P1::AbstractHPolygon{N},
                 P2::AbstractHPolygon{N},
                 [prune]::Bool=true
                ) where {N<:Real}

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
function intersection(P1::AbstractHPolygon{N},
                      P2::AbstractHPolygon{N},
                      prune::Bool=true
                     ) where {N<:Real}
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

    P = HPolygon(c, sort_constraints=false)
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
                 backend=default_lp_solver(N)) where {N<:Real}

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
`default_polyhedra_backend(P, N)` or use `CDDLib.Library()` for the CDD library.

There are some shortcomings of the removal of constraints using the default
Polyhedra library; see e.g. #1038 and Polyhedra#146. It is safer to check for
emptiness of intersection before calling this function in those cases.

### Algorithm

This implementation unifies the constraints of the two sets obtained from the
`constraints_list` method.
"""
function intersection(P1::AbstractPolyhedron{N},
                      P2::AbstractPolyhedron{N};
                      backend=default_lp_solver(N)) where {N<:Real}

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

# disambiguation
function intersection(S::AbstractSingleton{N},
                      P::AbstractPolyhedron{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, P)
end
function intersection(P::AbstractPolyhedron{N},
                      S::AbstractSingleton{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{typeof(S), LazySet{N}}, S, P)
end
function intersection(X::Interval{N}, P::AbstractPolyhedron{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{Interval{N}, LazySet{N}}, X, P)
end
function intersection(P::AbstractPolyhedron{N}, X::Interval{N}
                     ) where {N<:Real}
    return invoke(intersection, Tuple{Interval{N}, LazySet{N}}, X, P)
end

"""
    intersection(P1::Union{VPolytope{N}, VPolygon{N}},
                 P2::Union{VPolytope{N}, VPolygon{N}};
                 [backend]=default_polyhedra_backend(P1, N),
                 [prunefunc]=removevredundancy!) where {N<:Real}

Compute the intersection of two polytopes in vertex representation.

### Input

- `P1`        -- polytope in vertex representation
- `P2`        -- polytope in vertex representation
- `backend`   -- (optional, default: `default_polyhedra_backend(P1, N)`) the
                 backend for polyhedral computations
- `prunefunc` -- (optional, default: `removevredundancy!`) function to prune
                 the vertices of the result

### Output

A `VPolygon` if both arguments are `VPolygon`s, and a `VPolytope` otherwise.
"""
function intersection(P1::Union{VPolytope{N}, VPolygon{N}},
                      P2::Union{VPolytope{N}, VPolygon{N}};
                      backend=default_polyhedra_backend(P1, N),
                      prunefunc=removevredundancy!) where {N<:Real}
    Q1 = polyhedron(convert(VPolytope, P1); backend=backend)
    Q2 = polyhedron(convert(VPolytope, P2); backend=backend)
    Pint = Polyhedra.intersect(Q1, Q2)
    prunefunc(Pint)
    res = VPolytope(Pint)
    if P1 isa VPolygon && P2 isa VPolygon
        return convert(VPolygon, res)
    end
    return res
end

"""
    intersection(cup::UnionSet{N}, X::LazySet{N}) where {N<:Real}

Return the intersection of a union of two convex sets and another convex set.

### Input

- `cup` -- union of two convex sets
- `X`   -- convex set

### Output

The union of the pairwise intersections, expressed as a `UnionSet`.
If one of those sets is empty, only the other set is returned.
"""
function intersection(cup::UnionSet{N}, X::LazySet{N}) where {N<:Real}
    return intersection(cup.X, X) ∪ intersection(cup.Y, X)
end

# symmetric method
function intersection(X::LazySet{N}, cup::UnionSet{N}) where {N<:Real}
    return intersection(cup, X)
end

# disambiguation
function intersection(cup::UnionSet{N}, S::AbstractSingleton{N}) where {N<:Real}
    return element(S) ∈ cup ? S : EmptySet{N}(dim(S))
end
function intersection(S::AbstractSingleton{N}, cup::UnionSet{N}) where {N<:Real}
    return invoke(intersection, Tuple{UnionSet{N}, typeof(S)}, cup, S)
end

"""
    intersection(cup::UnionSetArray{N}, X::LazySet{N}) where {N<:Real}

Return the intersection of a union of a finite number of convex sets and another
convex set.

### Input

- `cup` -- union of a finite number of convex sets
- `X`   -- convex set

### Output

The union of the pairwise intersections, expressed as a `UnionSetArray`.
"""
function intersection(cup::UnionSetArray{N}, X::LazySet{N}) where {N<:Real}
    return UnionSetArray([intersection(Y, X) for Y in array(cup)])
end

# symmetric method
function intersection(X::LazySet{N}, cup::UnionSetArray{N}) where {N<:Real}
    return intersection(cup, X)
end

# disambiguation
function intersection(cup::UnionSetArray{N},
                      S::AbstractSingleton{N}) where {N<:Real}
    return element(S) ∈ cup ? S : EmptySet{N}(dim(S))
end
function intersection(S::AbstractSingleton{N},
                      cup::UnionSetArray{N}) where {N<:Real}
    return invoke(intersection, Tuple{UnionSetArray{N}, typeof(S)}, cup, S)
end

"""
    intersection(L::LinearMap{N}, S::LazySet{N}) where {N<:Real}

Return the intersection of a lazy linear map and a convex set.

### Input

- `L` -- linear map
- `S` -- convex set

### Output

The polytope obtained by the intersection of `l.M * L.X` and `S`.
"""
function intersection(L::LinearMap{N}, S::LazySet{N}) where {N<:Real}
    return intersection(linear_map(L.M, L.X), S)
end

# symmetric method
function intersection(S::LazySet{N}, L::LinearMap{N}) where {N}
    return intersection(L, S)
end

# disambiguation
function intersection(L1::LinearMap{N}, L2::LinearMap{N}) where {N}
    return intersection(linear_map(L1.M, L1.X), linear_map(L2.M, L2.X))
end

"""
    intersection(U::Universe{N}, X::LazySet{N}) where {N<:Real}

Return the intersection of a universe and a convex set.

### Input

- `U` -- universe
- `X` -- convex set

### Output

The set `X`.
"""
function intersection(U::Universe{N}, X::LazySet{N}) where {N<:Real}
    return X
end

# symmetric method
function intersection(X::LazySet{N}, U::Universe{N}) where {N<:Real}
    return X
end

# disambiguation
function intersection(U::Universe{N}, ::Universe{N}) where {N<:Real}
    return U
end
function intersection(U::Universe{N}, P::AbstractPolyhedron{N}) where {N<:Real}
    return P
end
function intersection(P::AbstractPolyhedron{N}, U::Universe{N}) where {N<:Real}
    return P
end
function intersection(U::Universe{N}, S::AbstractSingleton{N}) where {N<:Real}
    return S
end
function intersection(S::AbstractSingleton{N}, U::Universe{N}) where {N<:Real}
    return S
end
function intersection(X::Interval{N}, U::Universe{N}) where {N<:Real}
    return X
end
function intersection(U::Universe{N}, X::Interval{N}) where {N<:Real}
    return X
end

"""
    intersection(P::AbstractPolyhedron{N}, rm::ResetMap{N}) where {N<:Real}

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
function intersection(P::AbstractPolyhedron{N}, rm::ResetMap{N}) where {N<:Real}
    return intersection(P, HPolyhedron(constraints_list(rm)))
end

# symmetric method
function intersection(rm::ResetMap{N}, P::AbstractPolyhedron{N}) where {N<:Real}
    return intersection(P, rm)
end

# more efficient version for polytopic
function intersection(P::AbstractPolyhedron{N},
                      rm::ResetMap{N, <:AbstractPolytope}) where {N<:Real}
    return intersection(P, HPolytope(constraints_list(rm)))
end

# symmetric method
function intersection(rm::ResetMap{N, <:AbstractPolytope},
                      P::AbstractPolyhedron{N}) where {N<:Real}
    return intersection(P, rm)
end

function intersection(U::Universe{N}, X::CartesianProductArray{N}) where {N<:Real}
    return X
end

# symmetric method
function intersection(X::CartesianProductArray{N}, U::Universe{N}) where {N<:Real}
    return intersection(U, X)
end

"""
        intersection(X::CartesianProductArray{N}, Y::CartesianProductArray{N})

Return the intersection between cartesian products of a finite number of convex sets.

### Input

 - `X` -- cartesian product of a finite number of convex sets
 - `Y` -- cartesian product of a finite number of convex sets

### Output

The decomposed set which represents concrete intersection between `X` and `Y`

### Algorithm

This algorithm intersect corresponding blocks between sets.
"""
function intersection(X::CartesianProductArray{N}, Y::CartesianProductArray{N}) where {N<:Real}
    @assert same_block_structure(array(X), array(Y)) "block structure has to be the same"

    return CartesianProductArray([intersection(array(X)[i], array(Y)[i]) for i in eachindex(array(X))])
end

"""
    intersection(cpa::CartesianProductArray{N}, P::AbstractPolyhedron{N}
                ) where {N<:Real}

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
function intersection(cpa::CartesianProductArray{N}, P::AbstractPolyhedron{N}
                     ) where {N<:Real}
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
    cap = _intersection_cpa_polyhedron(
        CartesianProductArray(blocks[cb_start:cb_end]), P, dim_start:dim_end)

    # construct result
    result_array = vcat(blocks[1:(cb_start - 1)],  # unconstrained
                        [cap],  # constrained, intersected
                        blocks[(cb_end + 1):length(blocks)])  # unconstrained
    return CartesianProductArray(result_array)
end

function _intersection_cpa_polyhedron(cpa::CartesianProductArray{N},
                                      P::AbstractPolyhedron{N},
                                      constrained_dims) where {N<:Real}
    T = isbounded(cpa) ? HPolytope : HPolyhedron
    hpoly_low_dim = T(constraints_list(cpa))
    cap_low_dim = intersection(hpoly_low_dim, project(P, constrained_dims))
end

# symmetric method
function intersection(P::AbstractPolyhedron{N}, cpa::CartesianProductArray{N}
                     ) where {N<:Real}
    return intersection(cpa, P)
end
