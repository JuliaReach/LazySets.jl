import Base.issubset

"""
    issubset(X::LazySet, Y::LazySet, [witness]::Bool=false, args...)

Alias for `⊆` (inclusion check).

### Input

- `X`       -- set
- `Y`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ⊆ Y``
* If `witness` option is activated:
  * `(true, [])` iff ``X ⊆ Y``
  * `(false, v)` iff ``X ⊈ Y`` and ``v ∈ X \\setminus Y``

### Notes

For more documentation see `⊆`.
"""
function issubset end

# this operation is invalid, but it is a common error, so we give a detailed
# error message
function ⊆(::AbstractVector, ::LazySet)
    throw(ArgumentError("cannot make an inclusion check if the left-hand side " *
          "is a vector; either wrap it as a set with one element, as in " *
          "`Singleton(v) ⊆ X`, or check for set membership, as in `v ∈ X` " *
          "(they behave equivalently although the implementations may differ)"))
end

# conversion for IA types
⊆(X::LazySet, Y::IA.Interval) = ⊆(X, Interval(Y))
⊆(X::IA.Interval, Y::LazySet) = ⊆(Interval(X), Y)

⊆(X::LazySet, Y::IA.IntervalBox) = ⊆(X, convert(Hyperrectangle, Y))
⊆(X::IA.IntervalBox, Y::LazySet) = ⊆(convert(Hyperrectangle, X), Y)

"""
    ⊆(X::LazySet, P::LazySet, [witness]::Bool=false)

Check whether a set is contained in a polyhedral set, and if not, optionally
compute a witness.

### Input

- `X`       -- inner set
- `Y`       -- outer polyhedral set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ⊆ P``
* If `witness` option is activated:
  * `(true, [])` iff ``X ⊆ P``
  * `(false, v)` iff ``X ⊈ P`` and ``v ∈ X \\setminus P``

### Notes

We require that `constraints_list(P)` is available.

### Algorithm

We check inclusion of `X` in every constraint of `P`.
"""
function ⊆(X::LazySet, P::LazySet, witness::Bool=false)
    if is_polyhedral(P)
        return _issubset_constraints_list(X, P, witness)
    else
        error("an inclusion check for the given combination of set types is " *
              "not available")
    end
end

"""
    ⊆(S::LazySet, H::AbstractHyperrectangle, [witness]::Bool=false)

Check whether a set is contained in a hyperrectangular set, and if not,
optionally compute a witness.

### Input

- `S` -- inner set
- `H` -- outer hyperrectangular set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``S ⊆ H``
  * `(false, v)` iff ``S ⊈ H`` and ``v ∈ S \\setminus H``

### Algorithm

``S ⊆ H`` iff ``\\operatorname{ihull}(S) ⊆ H``, where  ``\\operatorname{ihull}``
is the interval-hull operator.
"""
function ⊆(S::LazySet, H::AbstractHyperrectangle, witness::Bool=false)
    return _issubset_in_hyperrectangle(S, H, witness)
end

function _issubset_in_hyperrectangle(S, H, witness)
    n = dim(S)
    @assert n == dim(H)
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
    if witness
        return (true, N[])
    else
        return true
    end
end

# disambiguation
⊆(P::AbstractPolytope, H::AbstractHyperrectangle, witness::Bool=false) =
    _issubset_in_hyperrectangle(P, H, witness)

"""
    ⊆(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle,
      [witness]::Bool=false)

Check whether a given hyperrectangular set is contained in another
hyperrectangular set, and if not, optionally compute a witness.

### Input

- `H1` -- inner hyperrectangular set
- `H2` -- outer hyperrectangular set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``H1 ⊆ H2``
* If `witness` option is activated:
  * `(true, [])` iff ``H1 ⊆ H2``
  * `(false, v)` iff ``H1 ⊈ H2`` and ``v ∈ H1 \\setminus H2``

### Algorithm

``H1 ⊆ H2`` iff ``c_1 + r_1 ≤ c_2 + r_2 ∧ c_1 - r_1 ≥ c_2 - r_2`` iff
``r_1 - r_2 ≤ c_1 - c_2 ≤ -(r_1 - r_2)``, where ``≤`` is taken component-wise.
"""
function ⊆(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle,
           witness::Bool=false)
    @assert dim(H1) == dim(H2)
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

    if witness
        return (true, N[])
    else
        return true
    end
end

"""
    ⊆(P::AbstractPolytope, S::LazySet, [witness]::Bool=false;
      [algorithm]="constraints")

Check whether a polytopic set is contained in a convex set, and if not,
optionally compute a witness.

### Input

- `P` -- inner polytopic set
- `S` -- outer convex set
- `witness`   -- (optional, default: `false`) compute a witness if activated
- `algorithm` -- (optional, default: `"constraints"`) algorithm for the
                 inclusion check; available options are:

    * `"constraints"`, using the list of constraints of `S` (requires that `S`
      is polyhedral) and support-function evaluations of `S`

    * `"vertices"`, using the list of vertices of `P` and membership evaluations
      of `S`

### Output

* If `witness` option is deactivated: `true` iff ``P ⊆ S``
* If `witness` option is activated:
  * `(true, [])` iff ``P ⊆ S``
  * `(false, v)` iff ``P ⊈ S`` and ``v ∈ P \\setminus S``

### Algorithm

- `"vertices"`:
Since ``S`` is convex, ``P ⊆ S`` iff ``v ∈ S`` for all vertices ``v`` of ``P``.
"""
function ⊆(P::AbstractPolytope, S::LazySet, witness::Bool=false;
           algorithm=nothing)
    @assert dim(P) == dim(S)
    if !isconvextype(typeof(S))
        error("an inclusion check for the given combination of set types is " *
              "not available")
    end

    if isnothing(algorithm)
        # TODO smarter evaluation which representation is better
        if is_polyhedral(S)
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
            if witness
                return (false, v)
            else
                return false
            end
        end
    end
    if witness
        return (true, N[])
    else
        return true
    end
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
  * `(false, v)` iff ``X ⊈ P`` and ``v ∈ P \\setminus X``

### Algorithm

Since ``X`` is convex, we can compare the support function of ``X`` and ``P`` in
each direction of the constraints of ``P``.

For witness generation, we use a support vector in the first direction where the
above check fails.
"""
function ⊆(X::LazySet, P::AbstractPolyhedron, witness::Bool=false)
    if !isconvextype(typeof(X))
        error("an inclusion check for the given combination of set types is " *
              "not available")
    end
    return _issubset_constraints_list(X, P, witness)
end

# for documentation see
# ⊆(X::LazySet, P::AbstractPolyhedron, witness::Bool=false)
function _issubset_constraints_list(S::LazySet, P::LazySet, witness::Bool=false)
    @assert dim(S) == dim(P)

    @inbounds for H in constraints_list(P)
        if !_leq(ρ(H.a, S), H.b)
            if witness
                return (false, σ(H.a, S))
            else
                return false
            end
        end
    end
    N = promote_type(eltype(S), eltype(P))
    return witness ? (true, N[]) : true
end

# disambiguations
for ST in [:AbstractPolytope, :AbstractHyperrectangle, :LineSegment]
    @eval ⊆(X::($ST), P::AbstractPolyhedron, witness::Bool=false) =
        _issubset_constraints_list(X, P, witness)
end

"""
    ⊆(S::AbstractSingleton, X::LazySet, [witness]::Bool=false)

Check whether a given set with a single value is contained in another set, and
if not, optionally compute a witness.

### Input

- `S`       -- inner set with a single value
- `X`       -- outer set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ⊆ X``
* If `witness` option is activated:
  * `(true, [])` iff ``S ⊆ X``
  * `(false, v)` iff ``S ⊈ X`` and ``v ∈ S \\setminus X``
"""
function ⊆(S::AbstractSingleton, X::LazySet, witness::Bool=false)
    return _issubset_singleton(S, X, witness)
end

function _issubset_singleton(S, X, witness)
    result = element(S) ∈ X

    if witness
        N = promote_type(eltype(S), eltype(X))
        return (result, result ? N[] : element(S))
    else
        return result
    end
end

# disambiguations
for ST in [:AbstractHyperrectangle, :AbstractPolyhedron, :UnionSetArray,
           :Complement]
    @eval ⊆(X::AbstractSingleton, Y::($ST), witness::Bool=false) =
        _issubset_singleton(X, Y, witness)
end

"""
    ⊆(S1::AbstractSingleton, S2::AbstractSingleton, witness::Bool=false)

Check whether a given set with a single value is contained in another set with a
single value, and if not, optionally compute a witness.

### Input

- `S1` -- inner set with a single value
- `S2` -- outer set with a single value
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S1 ⊆ S2`` iff ``S1 == S2``
* If `witness` option is activated:
  * `(true, [])` iff ``S1 ⊆ S2``
  * `(false, v)` iff ``S1 ⊈ S2`` and ``v ∈ S1 \\setminus S2``
"""
function ⊆(S1::AbstractSingleton, S2::AbstractSingleton, witness::Bool=false)
    result = _isapprox(element(S1), element(S2))
    if witness
        N = promote_type(eltype(S1), eltype(S2))
        return (result, result ? N[] : element(S1))
    else
        return result
    end
end

"""
    ⊆(B1::Ball2, B2::Ball2, [witness]::Bool=false)

Check whether a ball in the 2-norm is contained in another ball in the 2-norm,
and if not, optionally compute a witness.

### Input

- `B1` -- inner ball in the 2-norm
- `B2` -- outer ball in the 2-norm
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``B1 ⊆ B2``
* If `witness` option is activated:
  * `(true, [])` iff ``B1 ⊆ B2``
  * `(false, v)` iff ``B1 ⊈ B2`` and ``v ∈ B1 \\setminus B2``

### Algorithm

``B1 ⊆ B2`` iff ``‖ c_1 - c_2 ‖_2 + r_1 ≤ r_2``
"""
function ⊆(B1::Ball2, B2::Ball2, witness::Bool=false)
    result = norm(B1.center - B2.center, 2) + B1.radius <= B2.radius
    if !witness
        return result
    end

    if result
        N = promote_type(eltype(B1), eltype(B2))
        return (true, N[])
    end

    # compute a witness v
    v = B1.center .+ B1.radius * (B1.center .- B2.center)
    return (false, v)
end

"""
    ⊆(B::Union{Ball2, Ballp}, S::AbstractSingleton, witness::Bool=false)

Check whether a ball in the 2-norm or p-norm is contained in a set with a single
value, and if not, optionally compute a witness.

### Input

- `B` -- inner ball in the 2-norm or p-norm
- `S` -- outer set with a single value
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``B ⊆ S``
* If `witness` option is activated:
  * `(true, [])` iff ``B ⊆ S``
  * `(false, v)` iff ``B ⊈ S`` and ``v ∈ B \\setminus S``
"""
function ⊆(B::Union{Ball2, Ballp}, S::AbstractSingleton, witness::Bool=false)
    result = isapproxzero(B.radius) && _isapprox(B.center, element(S))
    if !witness
        return result
    end

    if result
        N = promote_type(eltype(B), eltype(S))
        return (result, N[])
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
    ⊆(L::LineSegment, S::LazySet, witness::Bool=false)

Check whether a line segment is contained in a convex set, and if not,
optionally compute a witness.

### Input

- `L` -- inner line segment
- `S` -- outer convex set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``L ⊆ S``
* If `witness` option is activated:
  * `(true, [])` iff ``L ⊆ S``
  * `(false, v)` iff ``L ⊈ S`` and ``v ∈ L \\setminus S``

### Algorithm

Since ``S`` is convex, ``L ⊆ S`` iff ``p ∈ S`` and ``q ∈ S``, where ``p, q`` are
the end points of ``L``.
"""
function ⊆(L::LineSegment, S::LazySet, witness::Bool=false)
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
    if !witness
        return result
    elseif result
        N = promote_type(eltype(L), eltype(S))
        return (result, N[])
    else
        return (result, p_in_S ? L.q : L.p)
    end
end

# disambiguation
⊆(L::LineSegment, H::AbstractHyperrectangle, witness::Bool=false) =
    _issubset_line_segment(L, H, witness)

"""
    ⊆(x::Interval, y::Interval, [witness]::Bool=false)

Check whether an interval is contained in another interval.

### Input

- `x`       -- inner interval
- `y`       -- outer interval
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

`true` iff ``x ⊆ y``.
"""
function ⊆(x::Interval, y::Interval, witness::Bool=false)
    if min(y) > min(x)
        witness && return (false, low(x))
        return false
    elseif max(x) > max(y)
        witness && return (false, high(x))
        return false
    end
    if !witness
        return true
    end
    N = promote_type(eltype(x), eltype(y))
    return (true, N[])
end

"""
    ⊆(x::Interval, U::UnionSet, [witness]::Bool=false)

Check whether an interval is contained in the union of two convex sets.

### Input

- `x` -- inner interval
- `U` -- outer union of two convex sets

### Output

`true` iff `x ⊆ U`.

### Notes

This implementation assumes that `U` is a union of one-dimensional convex sets.
Since these are equivalent to intervals, we convert to `Interval`s.

### Algorithm

Let ``U = a ∪ b `` where ``a`` and ``b`` are intervals and assume that the lower
bound of ``a`` is to the left of ``b``.
If the lower bound of ``x`` is to the left of ``a``, we have a counterexample.
Otherwise we compute the set difference ``y = x \\ a`` and check whether
``y ⊆ b`` holds.
"""
function ⊆(x::Interval, U::UnionSet, witness::Bool=false)
    @assert dim(U) == 1 "an interval is incompatible with a set of dimension " *
        "$(dim(U))"
    if !isconvextype(typeof(U.X)) || !isconvextype(typeof(U.Y))
        error("an inclusion check for the given combination of set types is " *
              "not available")
    end
    return _issubset_interval(x, convert(Interval, U.X), convert(Interval, U.Y),
                              witness)
end

function ⊆(x::Interval, U::UnionSet{N, <:Interval, <:Interval},
           witness::Bool=false) where {N}
    return _issubset_interval(x, U.X, U.Y, witness)
end

function _issubset_interval(x::Interval{N}, a::Interval, b::Interval,
                            witness) where {N}
    if min(a) > min(b)
        c = b
        b = a
        a = c
    end
    # a is on the left of b
    if min(x) < min(a)
        witness && return (false, low(x))
        return false
    end
    y = difference(x, a)
    if y ⊆ b
        witness && return (true, N[])
        return true
    elseif witness
        if min(b) > min(y)
            w = [(min(y) + min(b)) / 2]
        else
            w = high(y)
        end
        return (false, w)
    else
        return false
    end
end

function ⊆(x::Interval, U::UnionSetArray, witness::Bool=false)
    @assert dim(U) == 1 "an interval is incompatible with a set of dimension " *
        "$(dim(U))"
    V = _get_interval_array_copy(U)
    return _issubset_interval!(x, V, witness)
end

# disambiguation
function ⊆(x::Interval, U::UnionSetArray{N, <:AbstractHyperrectangle},
           witness::Bool=false) where {N}
    @assert dim(U) == 1 "an interval is incompatible with a set of dimension " *
        "$(dim(U))"
    V = _get_interval_array_copy(U)
    return _issubset_interval!(x, V, witness)
end

function _get_interval_array_copy(U::UnionSetArray)
    return [convert(Interval, X) for X in array(U)]
end

function _get_interval_array_copy(U::UnionSetArray{N, <:AbstractVector{<:Interval}}) where {N}
    return copy(array(U))
end

function _issubset_interval!(x::Interval{N}, intervals, witness) where {N}
    # sort intervals by lower bound
    sort!(intervals, lt=(x, y) -> min(x) <= min(y))

    # subtract intervals from x
    for y in intervals
        if min(y) > min(x)
            # lowest point of x is not contained
            witness && return (false, center(Interval(min(x), min(y))))
            return false
        end
        x = difference(x, y)
        if isempty(x)
            witness && return (true, N[])
            return true
        end
    end

    witness && return (false, low(x))
    return false
end

"""
    ⊆(X::LazySet{N}, U::UnionSetArray, witness::Bool=false;
      filter_redundant_sets::Bool=true) where {N}

Check whether a set is contained in a union of a finite number of sets.

### Input

- `X`       -- inner set
- `U`       -- outer union of a finite number of sets
- `witness` -- (optional, default: `false`) compute a witness if activated
- `filter_redundant_sets` -- (optional, default: `true`) ignore sets in `U` that
               do not intersect with `X`

### Output

`true` iff ``X ⊆ U``.

### Algorithm

This implementation is general and successively removes parts from `X` that are
covered by the sets in the union ``U`` using the `difference` function. For the
resulting subsets, this implementation relies on the methods `isdisjoint`, `⊆`,
and `intersection`.

As a preprocessing, this implementation checks if `X` is contained in any of the
sets in `U`.

The `filter_redundant_sets` option controls whether sets in `U` that do not
intersect with `X` should be ignored.
"""
function ⊆(X::LazySet, U::UnionSetArray, witness::Bool=false;
           filter_redundant_sets::Bool=true)
    return _issubset_unionsetarray(X, U, witness;
                                   filter_redundant_sets=filter_redundant_sets)
end

function _issubset_unionsetarray(X, U, witness::Bool=false;
                                 filter_redundant_sets::Bool=true)
    N = promote_type(eltype(X), eltype(U))

    # heuristics (necessary check): is X contained in any set in U?
    for rhs in array(U)
        if X ⊆ rhs
            if witness
                return (true, N[])
            else
                return true
            end
        end
    end

    if filter_redundant_sets
        # filter out those sets in U that do not intersect with X
        sets = Vector{eltype(array(U))}()
        for rhs in array(U)
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
            if witness
                return (false, an_element(Y))
            else
                return false
            end
        end
        # split wrt. the i-th set
        rhs = sets[idx]
        append!(queue, array(difference(Y, rhs)))
    end

    if witness
        return (true, N[])
    else
        return true
    end
end

# general container
function _inclusion_in_union_container(X::LazySet{N}, U::UnionSetArray) where {N}
    return LazySet{N}[]
end

# hyperrectangle container
function _inclusion_in_union_container(H::AbstractHyperrectangle,
                                       U::UnionSetArray{N, <:AbstractHyperrectangle}) where {N}
    return AbstractHyperrectangle{N}[]
end

# generally ignore check for flat sets
function _inclusion_in_union_isflat(X)
    return false
end

function _inclusion_in_union_isflat(H::AbstractHyperrectangle)
    return isflat(H)
end

# disambiguations
for ST in [:LineSegment]
    @eval ⊆(X::($ST), U::UnionSetArray, witness::Bool=false) =
        _issubset_unionsetarray(X, U, witness)
end
# disambiguation with additional kwarg
⊆(X::AbstractPolytope, U::UnionSetArray, witness::Bool=false; algorithm=nothing) =
    _issubset_unionsetarray(X, U, witness)


"""
    ⊆(∅::EmptySet, X::LazySet, witness::Bool=false)

Check whether the empty set is contained in another set.

### Input

- `∅`       -- inner empty set
- `X`       -- outer set
- `witness` -- (optional, default: `false`) compute a witness if activated
               (ignored, just kept for interface reasons)

### Output

`true`.
"""
function ⊆(∅::EmptySet, X::LazySet, witness::Bool=false)
    return _issubset_emptyset(∅, X, witness)
end

function _issubset_emptyset(∅::EmptySet, X::LazySet, witness::Bool=false)
    N = promote_type(eltype(∅), eltype(X))
    return witness ? (true, N[]) : true
end

# disambiguations
for ST in [:AbstractPolyhedron, :AbstractHyperrectangle, :Complement, :UnionSet,
           :UnionSetArray]
    @eval ⊆(∅::EmptySet, X::($ST), witness::Bool=false) =
        _issubset_emptyset(∅, X, witness)
end

"""
    ⊆(X::LazySet, ∅::EmptySet, [witness]::Bool=false)

Check whether a set is contained in the empty set.

### Input

- `X`       -- inner set
- `∅`       -- outer empty set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

`true` iff `X` is empty.

### Algorithm

We rely on `isempty(X)` for the emptiness check and on `an_element(X)` for
witness production.
"""
function ⊆(X::LazySet, ∅::EmptySet, witness::Bool=false)
    return _issubset_in_emptyset(X, ∅, witness)
end

function _issubset_in_emptyset(X::LazySet, ∅::EmptySet, witness::Bool=false)
    if isempty(X)
        N = promote_type(eltype(∅), eltype(X))
        return witness ? (true, N[]) : true
    else
        return witness ? (false, an_element(X)) : false
    end
end

# disambiguations
for ST in [:AbstractPolytope, :UnionSet, :UnionSetArray]
    @eval ⊆(X::($ST), ∅::EmptySet, witness::Bool=false) =
        _issubset_in_emptyset(X, ∅, witness)
end

# disambiguations for sets that are never empty
for ST in [:AbstractSingleton, :LineSegment]
    @eval ⊆(X::($ST), ::EmptySet, witness::Bool=false) =
        witness ? (false, an_element(X)) : false
end

function ⊆(∅₁::EmptySet, ∅₂::EmptySet, witness::Bool=false)
    N = promote_type(eltype(∅₁), eltype(∅₂))
    return witness ? (true, N[]) : true
end

"""
    ⊆(U::UnionSet, X::LazySet, [witness]::Bool=false)

Check whether a union of two convex sets is contained in another set.

### Input

- `U`       -- inner union of two convex sets
- `X`       -- outer set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``\\text{U} ⊆ X``
* If `witness` option is activated:
  * `(true, [])` iff ``\\text{U} ⊆ X``
  * `(false, v)` iff ``\\text{U} \\not\\subseteq X`` and
    ``v ∈ \\text{U} \\setminus X``
"""
function ⊆(U::UnionSet, X::LazySet, witness::Bool=false)
    return _issubset_union_in_set((U.X, U.Y), X, witness)
end

"""
    ⊆(U::UnionSetArray, X::LazySet, [witness]::Bool=false)

Check whether a union of a finite number of convex sets is contained in another
set.

### Input

- `U`       -- inner union of a finite number of convex sets
- `X`       -- outer set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``\\text{U} ⊆ X``
* If `witness` option is activated:
  * `(true, [])` iff ``\\text{U} ⊆ X``
  * `(false, v)` iff ``\\text{U} \\not\\subseteq X`` and
    ``v ∈ \\text{U} \\setminus X``
"""
function ⊆(U::UnionSetArray, X::LazySet, witness::Bool=false)
    return _issubset_union_in_set(array(U), X, witness)
end

# check for each set in `sets` that they are included in X
function _issubset_union_in_set(sets, X::LazySet{N}, witness::Bool=false) where {N}
    result = true
    v = N[]
    for Y in sets
        if witness
            result, v = ⊆(Y, X, witness)
        else
            result = ⊆(Y, X, witness)
        end
        if !result
            break
        end
    end
    return witness ? (result, v) : result
end

# disambiguations
for ST in [:AbstractHyperrectangle, :AbstractPolyhedron, :UnionSet,
           :UnionSetArray]
    @eval ⊆(U::UnionSet, X::($ST), witness::Bool=false) =
        _issubset_union_in_set((U.X, U.Y), X, witness)

    @eval ⊆(U::UnionSetArray, X::($ST), witness::Bool=false) =
        _issubset_union_in_set(array(U), X, witness)
end

"""
    ⊆(X::LazySet, U::Universe, [witness]::Bool=false)

Check whether a set is contained in a universe.

### Input

- `X`       -- inner set
- `U`       -- outer universe
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true`
* If `witness` option is activated: `(true, [])`
"""
function ⊆(X::LazySet, U::Universe, witness::Bool=false)
    return _issubset_universe(X, U, witness)
end

function _issubset_universe(X::LazySet, U::Universe, witness::Bool=false)
    if !witness
        return true
    end
    N = promote_type(eltype(X), eltype(U))
    return (true, N[])
end

# disambiguations
for ST in [:AbstractPolytope, :AbstractHyperrectangle, :AbstractSingleton,
           :LineSegment, :EmptySet, :UnionSet, :UnionSetArray]
    @eval ⊆(X::($ST), U::Universe, witness::Bool=false) =
        _issubset_universe(X, U, witness)
end

"""
    ⊆(U::Universe, X::LazySet, [witness]::Bool=false)

Check whether a universe is contained in another set, and otherwise optionally
compute a witness.

### Input

- `U`       -- inner universe
- `X`       -- outer set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``U ⊆ X``
* If `witness` option is activated:
  * `(true, [])` iff ``U ⊆ X``
  * `(false, v)` iff ``U \\not\\subseteq X`` and
    ``v ∈ U \\setminus X``

### Algorithm

We fall back to `isuniversal(X)`.
"""
function ⊆(::Universe, X::LazySet, witness::Bool=false)
    return isuniversal(X, witness)
end

# disambiguation
function ⊆(U1::Universe, U2::Universe, witness::Bool=false)
    return _issubset_universe(U1, U2, witness)
end

# disambiguations
for ST in [:AbstractPolyhedron, :AbstractPolytope, :AbstractHyperrectangle,
           :AbstractSingleton, :EmptySet, :UnionSetArray, :Complement]
    @eval ⊆(::Universe, X::($ST), witness::Bool=false) = isuniversal(X, witness)
end

"""
    ⊆(X::LazySet, C::Complement, [witness]::Bool=false)

Check whether a set is contained in the complement of another set, and otherwise
optionally compute a witness.

### Input

- `X`       -- inner set
- `C`       -- outer complement of a set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ⊆ C``
* If `witness` option is activated:
  * `(true, [])` iff ``X ⊆ C``
  * `(false, v)` iff ``X \\not\\subseteq C`` and
    ``v ∈ X \\setminus C``

### Algorithm

We fall back to `isdisjoint(X, C.X)`, which can be justified as follows.

```math
    X ⊆ Y^C ⟺ X ∩ Y = ∅
```
"""
function ⊆(X::LazySet, C::Complement, witness::Bool=false)
    return isdisjoint(X, C.X, witness)
end

# disambiguations
for ST in [:AbstractPolytope, :LineSegment, :UnionSet, :UnionSetArray]
    @eval ⊆(X::($ST), C::Complement, witness::Bool=false) =
        isdisjoint(X, C.X, witness)
end

"""
    ⊆(X::CartesianProduct, Y::CartesianProduct, [witness]::Bool=false;
      check_block_equality::Bool=true)

Check whether a Cartesian product of two sets is contained in another Cartesian
product of two sets, and otherwise optionally compute a witness.

### Input

- `X`       -- inner Cartesian product of two sets
- `Y`       -- outer Cartesian product of two sets
- `witness` -- (optional, default: `false`) compute a witness if activated
- `check_block_equality` -- (optional, default: `true`) flag for checking that
                            the block structure of the two sets is identical

### Output

* If `witness` option is deactivated: `true` iff ``X ⊆ Y``
* If `witness` option is activated:
  * `(true, [])` iff ``X ⊆ Y``
  * `(false, v)` iff ``X \\not\\subseteq Y`` and
    ``v ∈ X \\setminus Y``

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
function ⊆(X::CartesianProduct, Y::CartesianProduct, witness::Bool=false;
           check_block_equality::Bool=true)
    n1 = dim(X.X)
    n2 = dim(X.Y)
    if check_block_equality && (n1 != dim(Y.X) || n2 != dim(Y.Y))
        return _issubset_constraints_list(X, Y, witness)
    end

    # check first block
    result = ⊆(X.X, Y.X, witness)
    if !witness && !result
        return false
    elseif witness && !result[1]
        # construct a witness
        w = vcat(result[2], an_element(X.Y))
        return (false, w)
    end

    # check second block
    result = ⊆(X.Y, Y.Y, witness)
    if !witness && !result
        return false
    elseif witness && !result[1]
        # construct a witness
        w = vcat(an_element(X.X), result[2])
        return (false, w)
    end

    if !witness
        return true
    end
    N = promote_type(eltype(X), eltype(Y))
    return (true, N[])
end

"""
    ⊆(X::CartesianProductArray, Y::CartesianProductArray, [witness]::Bool=false;
      check_block_equality::Bool=true)

Check whether a Cartesian product of finitely many sets is contained in another
Cartesian product of finitely many sets, and otherwise optionally compute a
witness.

### Input

- `X`       -- inner Cartesian product of finitely many sets
- `Y`       -- outer Cartesian product of finitely many sets
- `witness` -- (optional, default: `false`) compute a witness if activated
- `check_block_equality` -- (optional, default: `true`) flag for checking that
                             the block structure of the two sets is identical

### Output

* If `witness` option is deactivated: `true` iff ``X ⊆ Y``
* If `witness` option is activated:
  * `(true, [])` iff ``X ⊆ Y``
  * `(false, v)` iff ``X \\not\\subseteq Y`` and
    ``v ∈ X \\setminus Y``

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
function ⊆(X::CartesianProductArray, Y::CartesianProductArray,
           witness::Bool=false; check_block_equality::Bool=true)
    aX = array(X)
    aY = array(Y)
    if check_block_equality && !same_block_structure(aX, aY)
        return _issubset_constraints_list(X, Y, witness)
    end

    N = promote_type(eltype(X), eltype(Y))
    for i in 1:length(aX)
        result = ⊆(aX[i], aY[i], witness)
        if !witness && !result
            return false
        elseif witness && !result[1]
            # construct a witness
            w = Vector{N}(undef, dim(X))
            k = 1
            for j in 1:length(aX)
                Xj = aX[j]
                l = k + dim(Xj)
                w[k:l-1] = j == i ? result[2] : an_element(Xj)
                k = l
            end
            return (false, w)
        end
    end

    return witness ? (true, N[]) : true
end

"""
    ⊆(Z::AbstractZonotope, H::AbstractHyperrectangle, [witness]::Bool=false)

Check whether a zonotopic set is contained in a hyperrectangular set.

### Input

- `Z`       -- inner zonotopic set
- `H`       -- outer hyperrectangular set
- `witness` -- (optional, default: `false`) compute a witness if activated
               (currently not supported)

### Output

`true` iff ``Z ⊆ H``.

### Algorithm

The algorithm is based on Lemma 3.1 in [1].

[1] Mitchell, I. M., Budzis, J., & Bolyachevets, A. *Invariant, viability and
discriminating kernel under-approximation via zonotope scaling*. HSCC 2019.
"""
function ⊆(Z::AbstractZonotope, H::AbstractHyperrectangle, witness::Bool=false)
    witness && raise(ArgumentError("witness production is not supported yet"))

    c = center(Z)
    G = genmat(Z)
    n, m = size(G)
    N = promote_type(eltype(Z), eltype(H))
    @inbounds for i = 1:n
        aux = zero(N)
        for j in 1:m
            aux += abs(G[i, j])
        end
        ubound = c[i] + aux
        lbound = c[i] - aux
        if !_leq(ubound, high(H, i)) || !_geq(lbound, low(H, i))
            return false
        end
    end
    return true
end

for ST in (AbstractZonotope, AbstractSingleton, LineSegment)
    @eval function ⊆(Z::$(ST), C::CartesianProduct{N, <:LazySet, <:Universe}) where N
        X = C.X
        Zp = project(Z, 1:dim(X))
        return ⊆(Zp, X)
    end

    @eval function ⊆(Z::$(ST), C::CartesianProduct{N, <:Universe, <:LazySet}) where N
        Y = C.Y
        Zp = project(Z, dim(C.X)+1:dim(C))
        return ⊆(Zp, Y)
    end

    # disambiguation
    @eval function ⊆(Z::$(ST), C::CartesianProduct{N, <:Universe, <:Universe}) where N
        return true
    end
end
