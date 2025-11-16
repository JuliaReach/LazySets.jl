# this operation is invalid, but it is a common error, so we give a detailed
# error message
function ⊆(::AbstractVector, ::LazySet)
    throw(ArgumentError("cannot make an inclusion check if the left-hand side " *
                        "is a vector; either wrap it as a set with one element, as in " *
                        "`Singleton(v) ⊆ X`, or check for set membership, as in `v ∈ X` " *
                        "(they behave equivalently although the implementations may differ)"))
end

"""
# Extended help

    ⊆(X::LazySet, Y::LazySet, witness::Bool=false)

### Algorithm

The default implementation assumes that `Y` is polyhedral, i.e., that
`constraints_list(Y)` is available, and checks inclusion of `X` in every
constraint of `Y`.
"""
@validate function ⊆(X::LazySet, Y::LazySet, witness::Bool=false)
    if ispolyhedral(Y)
        return _issubset_constraints_list(X, Y, witness)
    else
        throw(ArgumentError("an inclusion check for the given combination of " *
                            "set types is not available"))
    end
end

"""
# Extended help

    ⊆(S::LazySet, H::AbstractHyperrectangle, [witness]::Bool=false)

### Algorithm

``S ⊆ H`` iff ``\\operatorname{ihull}(S) ⊆ H``, where  ``\\operatorname{ihull}``
is the interval-hull operator.
"""
@validate function ⊆(S::LazySet, H::AbstractHyperrectangle, witness::Bool=false)
    return _issubset_in_hyperrectangle(S, H, witness)
end

function _issubset_in_hyperrectangle(S, H, witness)
    n = dim(S)
    N = promote_type(eltype(S), eltype(H))
    for i in 1:n
        lS, hS = extrema(S, i)
        lH, hH = extrema(H, i)
        if !witness && (lS < lH || hS > hH)
            return false
        elseif lS < lH
            # outside in negative direction
            v = σ(SingleEntryVector(i, n, -one(N)), S)
            return (false, v)
        elseif hS > hH
            # outside in positive direction
            v = σ(SingleEntryVector(i, n, one(N)), S)
            return (false, v)
        end
    end
    return _witness_result_empty(witness, true, N)
end

"""
# Extended help

    ⊆(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle,
      [witness]::Bool=false)

### Algorithm

``H1 ⊆ H2`` iff ``c_1 + r_1 ≤ c_2 + r_2 ∧ c_1 - r_1 ≥ c_2 - r_2`` iff
``r_1 - r_2 ≤ c_1 - c_2 ≤ -(r_1 - r_2)``, where ``≤`` is taken component-wise.
"""
@validate function ⊆(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle,
                     witness::Bool=false)
    N = promote_type(eltype(H1), eltype(H2))
    @inbounds for i in 1:dim(H1)
        c_dist = center(H1, i) - center(H2, i)
        r_dist = radius_hyperrectangle(H1, i) - radius_hyperrectangle(H2, i)
        # check if c_dist is not in the interval [r_dist, -r_dist]
        if !_leq(r_dist, c_dist) || !_leq(c_dist, -r_dist)
            if witness
                # compute a witness v
                v = copy(center(H1))
                if c_dist >= zero(N)
                    v[i] += radius_hyperrectangle(H1, i)
                else
                    v[i] -= radius_hyperrectangle(H1, i)
                end
                return (false, v)
            else
                return false
            end
        end
    end
    return _witness_result_empty(witness, true, N)
end

"""
# Extended help

    ⊆(P::AbstractPolytope, S::LazySet, [witness]::Bool=false;
      [algorithm]="constraints")

### Input

- `algorithm` -- (optional, default: `"constraints"`) algorithm for the
                 inclusion check; available options are:

    * `"constraints"`, using the list of constraints of `S` (requires that `S`
      is polyhedral) and support-function evaluations of `S`

    * `"vertices"`, using the list of vertices of `P` and membership evaluations
      of `S`

### Notes

`S` is assumed to be convex, which is asserted via `isconvextype`.

### Algorithm

- `"vertices"`:
Since ``S`` is convex, ``P ⊆ S`` iff ``v ∈ S`` for all vertices ``v`` of ``P``.
"""
@validate function ⊆(P::AbstractPolytope, S::LazySet, witness::Bool=false;
                     algorithm=nothing)
    if !isconvextype(typeof(S))
        error("an inclusion check for the given combination of set types is " *
              "not available")
    end

    if isnothing(algorithm)
        # TODO smarter evaluation which representation is better
        if ispolyhedral(S)
            algorithm = "constraints"
        else
            algorithm = "vertices"
        end
    end

    if algorithm == "constraints"
        return _issubset_constraints_list(P, S, witness)
    elseif algorithm == "vertices"
        return _issubset_vertices_list(P, S, witness)
    else
        error("algorithm $algorithm unknown")
    end
end

# check whether P ⊆ S by testing whether each vertex of P belongs to S
function _issubset_vertices_list(P, S, witness)
    for v in vertices(P)
        if v ∉ S
            return _witness_result(witness, false, v)
        end
    end
    return _witness_result_empty(witness, true, P, S)
end

"""
    ⊆(X::LazySet, P::AbstractPolyhedron, [witness]::Bool=false)

Check whether a convex set is contained in a polyhedral set, and if not,
optionally compute a witness.

### Input

- `X` -- inner convex set
- `P` -- outer polyhedral set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ⊆ P``
* If `witness` option is activated:
  * `(true, [])` iff ``X ⊆ P``
  * `(false, v)` iff ``X ⊈ P`` and ``v ∈ P ∖ X``

### Algorithm

Since ``X`` is convex, we can compare the support function of ``X`` and ``P`` in
each direction of the constraints of ``P``.

For witness generation, we use a support vector in the first direction where the
above check fails.

See [WetzlingerKBA23; Proposition 7](@citet).
"""
@validate function ⊆(X::LazySet, P::AbstractPolyhedron, witness::Bool=false)
    if !isconvextype(typeof(X))
        return _issubset_in_polyhedron_high(X, P, witness)
    end
    return _issubset_constraints_list(X, P, witness)
end

# S ⊆ P where P = ⟨Cx ≤ d⟩  iff  y ≤ d where y is the upper corner of box(C*S).
# See [WetzlingerKBA23; Proposition 7](@citet).
function _issubset_in_polyhedron_high(S::LazySet, P::LazySet, witness::Bool=false)
    C, d = tosimplehrep(P)
    x = high(C * S)
    result = all(x .≤ d)

    if result
        return _witness_result_empty(witness, true, S, P)
    elseif !witness
        return false
    end
    throw(ArgumentError("witness production is not supported yet"))
end

@validate function ⊆(Z::AbstractZonotope, P::AbstractPolyhedron, witness::Bool=false)
    return _issubset_zonotope_in_polyhedron(Z, P, witness)
end

# See [WetzlingerKBA23; Proposition 7](@citet).
function _issubset_zonotope_in_polyhedron(Z::AbstractZonotope, P::LazySet,
                                          witness::Bool=false)
    # corner case: no generator
    if ngens(Z) == 0
        c = center(Z)
        result = c ∈ P
        return _witness_result_empty(witness, result, Z, P, c)
    end

    C, d = tosimplehrep(P)
    c = center(Z)
    G = genmat(Z)
    A = sum(abs.(C * gj) for gj in eachcol(G))
    b = d - C * c
    result = all(_leq.(A, b))

    if result
        return _witness_result_empty(witness, true, Z, P)
    elseif !witness
        return false
    end
    # fall back to default algorithm for witness production
    # TODO more efficient implementation
    return _issubset_constraints_list(Z, P, witness)
end

# for documentation see
# ⊆(X::LazySet, P::AbstractPolyhedron, witness::Bool=false)
function _issubset_constraints_list(S::LazySet, P::LazySet, witness::Bool=false)
    @assert ispolyhedral(P) "this inclusion check requires a polyhedral set " *
                            "on the right-hand side"

    @inbounds for H in constraints_list(P)
        if !_leq(ρ(H.a, S), H.b)
            return witness ? (false, σ(H.a, S)) : false
        end
    end
    return _witness_result_empty(witness, true, S, P)
end

@validate function ⊆(S1::AbstractSingleton, S2::AbstractSingleton, witness::Bool=false)
    s1 = element(S1)
    result = _isapprox(s1, element(S2))
    return _witness_result_empty(witness, result, S1, S2, s1)
end

@validate function ⊆(S::AbstractSingleton, X::LazySet, witness::Bool=false)
    return _issubset_singleton(S, X, witness)
end

function _issubset_singleton(S, X, witness)
    s = element(S)
    result = s ∈ X
    return _witness_result_empty(witness, result, S, X, s)
end

@validate function ⊆(B::Union{Ball2,Ballp}, S::AbstractSingleton, witness::Bool=false)
    result = isapproxzero(B.radius) && _isapprox(B.center, element(S))
    if result
        return _witness_result_empty(witness, true, B, S)
    elseif !witness
        return false
    end

    # compute a witness v
    if B.center != element(S)
        v = B.center
    else
        v = copy(B.center)
        v[1] += B.radius
    end
    return (false, v)
end

"""
# Extended help

    ⊆(L::LineSegment, S::LazySet, witness::Bool=false)

### Algorithm

Since ``S`` is convex, ``L ⊆ S`` iff ``p ∈ S`` and ``q ∈ S``, where ``p, q`` are
the end points of ``L``.
"""
@validate function ⊆(L::LineSegment, S::LazySet, witness::Bool=false)
    if !isconvextype(typeof(S))
        error("an inclusion check for the given combination of set types is " *
              "not available")
    end

    return _issubset_line_segment(L, S, witness)
end

# requires convexity of S
function _issubset_line_segment(L, S, witness)
    p_in_S = L.p ∈ S
    result = p_in_S && L.q ∈ S
    if result
        return _witness_result_empty(witness, true, L, S)
    end
    return witness ? (false, p_in_S ? L.q : L.p) : false
end

"""
    ⊆(X::Interval, U::UnionSet, [witness]::Bool=false)

Check whether an interval is contained in the union of two convex sets.

### Input

- `X` -- inner interval
- `U` -- outer union of two convex sets

### Output

`true` iff `X ⊆ U`.

### Notes

This implementation assumes that `U` is a union of one-dimensional convex sets.
Since these are equivalent to intervals, we convert to `Interval`s.

### Algorithm

Let ``U = Y ∪ Z `` where ``Y`` and ``Z`` are intervals and assume that the lower
bound of ``Y`` is to the left of ``Z``.
If the lower bound of ``X`` is to the left of ``Y``, we have a counterexample.
Otherwise we compute the set difference ``W = X \\ Y`` and check whether
``W ⊆ Z`` holds.
"""
@validate function ⊆(X::Interval, U::UnionSet, witness::Bool=false)
    if !isconvextype(typeof(first(U))) || !isconvextype(typeof(second(U)))
        error("an inclusion check for the given combination of set types is " *
              "not available")
    end
    return _issubset_interval(X, convert(Interval, first(U)),
                              convert(Interval, second(U)), witness)
end

@validate function ⊆(X::Interval, U::UnionSet{N,<:Interval,<:Interval},
                     witness::Bool=false) where {N}
    return _issubset_interval(X, first(U), second(U), witness)
end

function _issubset_interval(X::Interval{N}, Y::Interval, Z::Interval,
                            witness) where {N}
    if min(Y) > min(Z)
        W = Z
        Z = Y
        Y = W
    end
    # a is on the left of b
    if min(X) < min(Y)
        return witness ? (false, low(X)) : false
    end
    W = difference(X, Y)
    if W ⊆ Z
        return _witness_result_empty(witness, true, N)
    elseif !witness
        return false
    end

    # compute witness
    w = min(Z) > min(W) ? [(min(W) + min(Z)) / 2] : high(W)
    return (false, w)
end

@validate function ⊆(X::Interval, U::UnionSetArray, witness::Bool=false)
    V = _get_interval_array_copy(U)
    return _issubset_interval!(X, V, witness)
end

function _get_interval_array_copy(U::UnionSetArray{N}) where {N}
    out = Vector{LazySet{N}}(undef, length(array(U)))
    for (i, Xi) in enumerate(array(U))
        Yi = _to_unbounded_interval(Xi)
        if Yi isa Universe
            return Yi
        end
        out[i] = Yi
    end
    return out
end

_to_unbounded_interval(X::Interval) = X
_to_unbounded_interval(H::HalfSpace) = H
_to_unbounded_interval(U::Universe) = U

function _to_unbounded_interval(X::LazySet{N}) where {N}
    if !isconvextype(typeof(X))
        throw(ArgumentError("unions with non-convex sets are not supported"))
    end
    l, h = extrema(X, 1)
    if isinf(l)
        if isinf(h)
            return [Universe{N}(1)]
        else
            return HalfSpace([one(N)], h)
        end
    elseif isinf(h)
        return HalfSpace([-one(N)], -l)
    else
        return Interval(l, h)
    end
end

function _get_interval_array_copy(U::UnionSetArray{N,<:AbstractHyperrectangle}) where {N}
    return [convert(Interval, X) for X in U]
end

function _get_interval_array_copy(U::UnionSetArray{N,<:Interval}) where {N}
    return copy(array(U))
end

function _issubset_interval!(X::Interval{N}, intervals, witness) where {N}
    # sort intervals by lower bound
    sort!(intervals; lt=(x, y) -> low(x, 1) <= low(y, 1))

    # subtract intervals from x
    for Y in intervals
        if low(Y, 1) > low(X, 1)
            # lowest point of x is not contained
            # witness is the point in the middle
            return witness ? (false, [(low(X, 1) + low(Y, 1)) / 2]) : false
        end
        X = difference(X, Y)
        if isempty(X)
            return _witness_result_empty(witness, true, N)
        end
    end

    return witness ? (false, center(X)) : false
end

"""
# Extended help

    ⊆(X::LazySet{N}, U::UnionSetArray, witness::Bool=false;
      filter_redundant_sets::Bool=true) where {N}

### Input

- `filter_redundant_sets` -- (optional, default: `true`) ignore sets in `U` that
               do not intersect with `X`

### Algorithm

This implementation is general and successively removes parts from `X` that are
covered by the sets in the union ``U`` using the `difference` function. For the
resulting subsets, this implementation relies on the methods `isdisjoint`, `⊆`,
`intersection`, and `difference`.

As a preprocessing, this implementation checks if `X` is contained in any of the
sets in `U`.

The `filter_redundant_sets` option controls whether sets in `U` that do not
intersect with `X` should be ignored.
"""
@validate function ⊆(X::LazySet, U::UnionSetArray, witness::Bool=false;
                     filter_redundant_sets::Bool=true)
    return _issubset_unionsetarray(X, U, witness;
                                   filter_redundant_sets=filter_redundant_sets)
end

function _issubset_unionsetarray(X, U, witness::Bool=false;
                                 filter_redundant_sets::Bool=true)
    # heuristics (necessary check): is X contained in any set in U?
    for rhs in U
        if X ⊆ rhs
            return _witness_result_empty(witness, true, X, U)
        end
    end

    if filter_redundant_sets
        # filter out those sets in U that do not intersect with X
        sets = Vector{eltype(array(U))}()
        for rhs in U
            if !isdisjoint(X, rhs)
                push!(sets, rhs)
            end
        end
    else
        sets = array(U)
    end

    queue = _inclusion_in_union_container(X, U)
    push!(queue, X)
    while !isempty(queue)
        Y = pop!(queue)
        # first check if Y is fully contained to avoid splitting/recursion
        # keep track of the first set that intersects with Y
        idx = 0
        contained = false
        for (i, rhs) in enumerate(sets)
            if Y ⊆ rhs
                contained = true
                break
            elseif idx == 0 && !isdisjoint(Y, rhs)
                # does not terminate if the intersection is flat
                cap = intersection(Y, rhs)
                if !isempty(cap) && !_inclusion_in_union_isflat(cap)
                    idx = i
                end
            end
        end
        if contained
            continue
        end
        if idx == 0
            return witness ? (false, an_element(Y)) : false
        end
        # split wrt. the i-th set
        rhs = sets[idx]
        append!(queue, array(difference(Y, rhs)))
    end
    return _witness_result_empty(witness, true, X, U)
end

# general container
function _inclusion_in_union_container(X::LazySet{N}, U::UnionSetArray) where {N}
    return LazySet{N}[]
end

# hyperrectangle container
function _inclusion_in_union_container(H::AbstractHyperrectangle,
                                       U::UnionSetArray{N,<:AbstractHyperrectangle}) where {N}
    return AbstractHyperrectangle{N}[]
end

# generally ignore check for flat sets
function _inclusion_in_union_isflat(X)
    return false
end

function _inclusion_in_union_isflat(H::AbstractHyperrectangle)
    return isflat(H)
end

@validate function ⊆(∅::EmptySet, X::LazySet, witness::Bool=false)
    return _issubset_emptyset(∅, X, witness)
end

"""
# Extended help

    ⊆(X::LazySet, ∅::EmptySet, [witness]::Bool=false)

### Algorithm

We rely on `isempty(X)` for the emptiness check and on `an_element(X)` for
witness production.
"""
@validate function ⊆(X::LazySet, ∅::EmptySet, witness::Bool=false)
    return _issubset_emptyset2(X, ∅, witness)
end

"""
# Extended help

    ⊆(U::UnionSet, X::LazySet, [witness]::Bool=false)

### Notes

This implementation assumes that the sets in the union `U` are convex.
"""
@validate function ⊆(U::UnionSet, X::LazySet, witness::Bool=false)
    return _issubset_union_in_set(U, X, witness)
end

"""
# Extended help

    ⊆(U::UnionSetArray, X::LazySet, [witness]::Bool=false)

### Notes

This implementation assumes that the sets in the union `U` are convex.
"""
@validate function ⊆(U::UnionSetArray, X::LazySet, witness::Bool=false)
    return _issubset_union_in_set(U, X, witness)
end

# check for each set in `sets` that they are included in X
function _issubset_union_in_set(cup::Union{UnionSet,UnionSetArray}, X::LazySet{N},
                                witness::Bool=false) where {N}
    result = true
    v = N[]
    for Y in cup
        if witness
            result, v = ⊆(Y, X, witness)
        else
            result = ⊆(Y, X, witness)
        end
        if !result
            break
        end
    end
    return _witness_result(witness, result, v)
end

"""
# Extended help

    ⊆(U::Universe, X::LazySet, [witness]::Bool=false)

### Algorithm

We fall back to `isuniversal(X)`.
"""
@validate function ⊆(U::Universe, X::LazySet, witness::Bool=false)
    return _issubset_universe(U, X, witness)
end

@validate function ⊆(X::LazySet, U::Universe, witness::Bool=false)
    return _issubset_universe2(X, U, witness)
end

"""
# Extended help

    ⊆(X::LazySet, C::Complement, [witness]::Bool=false)

### Algorithm

We fall back to `isdisjoint(X, C.X)`, which can be justified as follows.

```math
    X ⊆ Y^C ⟺ X ∩ Y = ∅
```
"""
@validate function ⊆(X::LazySet, C::Complement, witness::Bool=false)
    return isdisjoint(X, C.X, witness)
end

"""
# Extended help

    ⊆(X::CartesianProduct, Y::CartesianProduct, [witness]::Bool=false;
      check_block_equality::Bool=true)

### Input

- `check_block_equality` -- (optional, default: `true`) flag for checking that
                            the block structure of the two sets is identical

### Notes

This algorithm requires that the two Cartesian products share the same block
structure.
If `check_block_equality` is activated, we check this property and, if it does
not hold, we use a fallback implementation based on conversion to constraint
representation (assuming that the sets are polyhedral).

### Algorithm

We check inclusion for each block of the Cartesian products.

For witness production, we obtain a witness in one of the blocks.
We then construct a high-dimensional witness by obtaining any point in the other
blocks (using `an_element`) and concatenating these (lower-dimensional) points.
"""
@validate function ⊆(X::CartesianProduct, Y::CartesianProduct, witness::Bool=false;
                     check_block_equality::Bool=true)
    n1 = dim(first(X))
    n2 = dim(second(X))
    if check_block_equality && (n1 != dim(first(Y)) || n2 != dim(second(Y)))
        return _issubset_constraints_list(X, Y, witness)
    end

    # check first block
    result = ⊆(first(X), first(Y), witness)
    if !witness && !result
        return false
    elseif witness && !result[1]
        # construct a witness
        w = vcat(result[2], an_element(second(X)))
        return (false, w)
    end

    # check second block
    result = ⊆(second(X), second(Y), witness)
    if !witness && !result
        return false
    elseif witness && !result[1]
        # construct a witness
        w = vcat(an_element(first(X)), result[2])
        return (false, w)
    end
    return _witness_result_empty(witness, true, X, Y)
end

"""
# Extended help

    ⊆(X::CartesianProductArray, Y::CartesianProductArray, [witness]::Bool=false;
      check_block_equality::Bool=true)

### Input

- `check_block_equality` -- (optional, default: `true`) flag for checking that
                             the block structure of the two sets is identical

### Notes

This algorithm requires that the two Cartesian products share the same block
structure.
If `check_block_equality` is activated, we check this property and, if it does
not hold, we use a fallback implementation based on conversion to constraint
representation (assuming that the sets are polyhedral).

### Algorithm

We check inclusion for each block of the Cartesian products.

For witness production, we obtain a witness in one of the blocks.
We then construct a high-dimensional witness by obtaining any point in the other
blocks (using `an_element`) and concatenating these (lower-dimensional) points.
"""
@validate function ⊆(X::CartesianProductArray, Y::CartesianProductArray,
                     witness::Bool=false; check_block_equality::Bool=true)
    aX = array(X)
    aY = array(Y)
    if check_block_equality && !same_block_structure(aX, aY)
        return _issubset_constraints_list(X, Y, witness)
    end

    N = promote_type(eltype(X), eltype(Y))
    for i in eachindex(aX)
        result = ⊆(aX[i], aY[i], witness)
        if !witness && !result
            return false
        elseif witness && !result[1]
            # construct a witness
            w = Vector{N}(undef, dim(X))
            k = 1
            for j in eachindex(aX)
                Xj = aX[j]
                l = k + dim(Xj)
                w[k:(l - 1)] = j == i ? result[2] : an_element(Xj)
                k = l
            end
            return (false, w)
        end
    end
    return _witness_result_empty(witness, true, N)
end

"""
# Extended help

    ⊆(Z::AbstractZonotope, H::AbstractHyperrectangle, [witness]::Bool=false)

### Algorithm

The algorithm is based on [MitchellBB19; Lemma 3.1](@citet).
"""
@validate function ⊆(Z::AbstractZonotope, H::AbstractHyperrectangle, witness::Bool=false)
    c = center(Z)
    G = genmat(Z)
    n, m = size(G)
    N = promote_type(eltype(Z), eltype(H))
    @inbounds for i in 1:n
        aux = zero(N)
        for j in 1:m
            aux += abs(G[i, j])
        end
        ubound = c[i] + aux
        lbound = c[i] - aux
        if !_leq(ubound, high(H, i)) || !_geq(lbound, low(H, i))
            if witness
                throw(ArgumentError("witness production is not supported yet"))
            end
            return false
        end
    end
    return _witness_result_empty(witness, true, Z, H)
end

for T in (AbstractZonotope, AbstractSingleton, LineSegment)
    @eval @validate function ⊆(Z::$(T), C::CartesianProduct{N,<:LazySet,<:Universe}) where {N}
        X = first(C)
        Zp = project(Z, 1:dim(X))
        return ⊆(Zp, X)
    end

    @eval @validate function ⊆(Z::$(T), C::CartesianProduct{N,<:Universe,<:LazySet}) where {N}
        Y = second(C)
        Zp = project(Z, (dim(first(C)) + 1):dim(C))
        return ⊆(Zp, Y)
    end

    # disambiguation
    @eval @validate function ⊆(X::$(T), Y::CartesianProduct{N,<:Universe,<:Universe}) where {N}
        return true
    end
end

# ============== #
# disambiguation #
# ============== #

@validate function ⊆(P::AbstractPolytope, H::AbstractHyperrectangle, witness::Bool=false)
    return _issubset_in_hyperrectangle(P, H, witness)
end

for T in (:AbstractPolytope, :AbstractHyperrectangle, :LineSegment)
    @eval @validate function ⊆(X::($T), P::AbstractPolyhedron, witness::Bool=false)
        return _issubset_constraints_list(X, P, witness)
    end
end

for T in (:AbstractHyperrectangle, :AbstractPolyhedron, :UnionSetArray, :Complement)
    @eval @validate function ⊆(X::AbstractSingleton, Y::($T), witness::Bool=false)
        return _issubset_singleton(X, Y, witness)
    end
end

@validate function ⊆(L::LineSegment, H::AbstractHyperrectangle, witness::Bool=false)
    return _issubset_line_segment(L, H, witness)
end

@validate function ⊆(X::Interval, U::UnionSetArray{N,<:AbstractHyperrectangle},
                     witness::Bool=false) where {N}
    V = _get_interval_array_copy(U)
    return _issubset_interval!(X, V, witness)
end

@validate function ⊆(X::LineSegment, U::UnionSetArray, witness::Bool=false)
    return _issubset_unionsetarray(X, U, witness)
end

# with additional kwarg
@validate function ⊆(X::AbstractPolytope, U::UnionSetArray, witness::Bool=false; algorithm=nothing)
    return _issubset_unionsetarray(X, U, witness)
end

for T in (:AbstractPolyhedron, :AbstractHyperrectangle, :Complement, :UnionSet, :UnionSetArray)
    @eval @validate function ⊆(∅::EmptySet, X::($T), witness::Bool=false)
        return _issubset_emptyset(∅, X, witness)
    end
end

for T in (:AbstractPolytope, :UnionSet, :UnionSetArray, :AbstractSingleton, :LineSegment)
    @eval @validate function ⊆(X::($T), ∅::EmptySet, witness::Bool=false)
        return _issubset_emptyset2(X, ∅, witness)
    end
end

for TU in (:UnionSet, :UnionSetArray)
    for T in (:AbstractHyperrectangle, :AbstractPolyhedron, :UnionSet, :UnionSetArray)
        @eval @validate function ⊆(U::($TU), X::($T), witness::Bool=false)
            return _issubset_union_in_set(U, X, witness)
        end
    end
end

for T in (:AbstractPolyhedron, :AbstractPolytope, :AbstractHyperrectangle,
          :AbstractSingleton, :EmptySet, :UnionSetArray, :Complement)
    @eval @validate function ⊆(U::Universe, X::($T), witness::Bool=false)
        return _issubset_universe(U, X, witness)
    end
end

for T in (:AbstractPolytope, :AbstractZonotope, :AbstractHyperrectangle,
          :AbstractSingleton, :LineSegment, :EmptySet, :UnionSet, :UnionSetArray)
    @eval @validate function ⊆(X::($T), U::Universe, witness::Bool=false)
        return _issubset_universe2(X, U, witness)
    end
end

for T in (:AbstractPolytope, :LineSegment, :UnionSet, :UnionSetArray)
    @eval @validate function ⊆(X::($T), C::Complement, witness::Bool=false)
        return isdisjoint(X, C.X, witness)
    end
end
