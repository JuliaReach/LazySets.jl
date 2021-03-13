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

# this operation is forbidden, but it is a common error so we give a detailed error message
function ⊆(::AbstractVector, ::LazySet)
    throw(ArgumentError("cannot make an inclusion check if the left-hand side " *
          "is a vector; either wrap it as a set with one element, as in " *
          "`Singleton(v) ⊆ X`, or check for set membership, as in `v ∈ X` " *
          "(they behave equivalently although the implementations may differ)"))
end


# --- fall-back with constraints_list of rhs ---

"""
    ⊆(X::LazySet, P::LazySet, [witness]::Bool=false)

Check whether a convex set is contained in a polyhedral set, and if not,
optionally compute a witness.

### Input

- `X`       -- inner convex set
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
    if applicable(constraints_list, P)
        return _issubset_constraints_list(X, P, witness)
    else
        error("an inclusion check for the given combination of set types is " *
              "not available")
    end
end


# --- AbstractHyperrectangle ---


"""
    ⊆(S::LazySet, H::AbstractHyperrectangle, [witness]::Bool=false)

Check whether a convex set is contained in a hyperrectangular set, and if not,
optionally compute a witness.

### Input

- `S` -- inner convex set
- `H` -- outer hyperrectangular set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``S ⊆ H``
  * `(false, v)` iff ``S ⊈ H`` and ``v ∈ S \\setminus H``

### Algorithm

``S ⊆ H`` iff ``\\operatorname{ihull}(S) ⊆ H``, where  ``\\operatorname{ihull}``
is the interval hull operator.
"""
function ⊆(S::LazySet, H::AbstractHyperrectangle, witness::Bool=false)
    return _issubset_in_hyperrectangle(S, H, witness)
end

function _issubset_in_hyperrectangle(S, H, witness)
    Shull = Approximations.interval_hull(S)
    return ⊆(Shull, H, witness)
end

# fix ambiguity
⊆(P::AbstractPolytope, H::AbstractHyperrectangle, witness::Bool=false) = _issubset_in_hyperrectangle(P, H, witness)

"""
    ⊆(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle, [witness]::Bool=false)

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
function ⊆(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle, witness::Bool=false)
    @assert dim(H1) == dim(H2)
    N = promote_type(eltype(H1), eltype(H2))

    @inbounds for i in 1:dim(H1)
        c_dist = center(H1, i) - center(H2, i)
        r_dist = radius_hyperrectangle(H1, i) - radius_hyperrectangle(H2, i)
        # check if c_dist is not in the interval [r_dist, -r_dist]
        if !_leq(r_dist, c_dist) || !_leq(c_dist, -r_dist)
            if witness
                # compute a witness 'p' in the difference
                p = copy(center(H1))
                if c_dist >= zero(N)
                    p[i] += radius_hyperrectangle(H1, i)
                else
                    p[i] -= radius_hyperrectangle(H1, i)
                end
                return (false, p)
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


# --- AbstractPolytope ---


"""
    ⊆(P::AbstractPolytope, S::LazySet, [witness]::Bool=false;
      algorithm=_default_issubset(P, S))

Check whether a polytope is contained in a convex set, and if not, optionally
compute a witness.

### Input

- `P` -- inner polytope
- `S` -- outer convex set
- `witness`   -- (optional, default: `false`) compute a witness if activated
- `algorithm` -- (optional, default: `"constraints"` if the constraints list of `S`
                 is available, otherwise `"vertices"`) algorithm for the inclusion
                 check; available options are:

    * `"constraints"`, using the list of constraints of `P` and support function
      evaluations of `S`

    * `"vertices"`, using the list of vertices of `P` and membership evaluations
      of `S`

### Output

* If `witness` option is deactivated: `true` iff ``P ⊆ S``
* If `witness` option is activated:
  * `(true, [])` iff ``P ⊆ S``
  * `(false, v)` iff ``P ⊈ S`` and ``v ∈ P \\setminus S``

### Algorithm

Since ``S`` is convex, ``P ⊆ S`` iff ``v_i ∈ S`` for all vertices ``v_i`` of
``P``.
"""
function ⊆(P::AbstractPolytope, S::LazySet, witness::Bool=false;
           algorithm=_default_issubset(P, S))
    @assert dim(P) == dim(S)

    if algorithm == "constraints"
        return _issubset_constraints_list(P, S, witness)
    elseif algorithm == "vertices"
        return _issubset_vertices_list(P, S, witness)
    else
        error("algorithm $algorithm unknown")
    end
end

@inline function _default_issubset(P, S)
    if applicable(constraints_list, S)
        return "constraints"
    else
        return "vertices"
    end
end

# check whether P ⊆ S by testing if each vertex of P belongs to S
function _issubset_vertices_list(P, S, witness)
    for v in vertices_list(P)
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

Check whether a convex set is contained in a polyhedron, and if not, optionally
compute a witness.

### Input

- `X` -- inner convex set
- `P` -- outer polyhedron (including a half-space)
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ⊆ P``
* If `witness` option is activated:
  * `(true, [])` iff ``X ⊆ P``
  * `(false, v)` iff ``X ⊈ P`` and ``v ∈ P \\setminus X``

### Algorithm

Since ``X`` is convex, we can compare the support function of ``X`` and ``P`` in
each direction of the constraints of ``P``.

For witness generation, we use the support vector in the first direction where
the above check fails.
"""
function ⊆(X::LazySet, P::AbstractPolyhedron, witness::Bool=false)
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
for ST in [AbstractPolytope, AbstractHyperrectangle, LineSegment]
    @eval ⊆(X::($ST), P::AbstractPolyhedron, witness::Bool=false) = _issubset_constraints_list(X, P, witness)
end

# --- AbstractSingleton ---

"""
    ⊆(S::AbstractSingleton, X::LazySet, [witness]::Bool=false)

Check whether a given set with a single value is contained in a convex set, and
if not, optionally compute a witness.

### Input

- `S`       -- inner set with a single value
- `X`       -- outer convex set
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
for ST in [AbstractHyperrectangle, AbstractPolyhedron]
    @eval ⊆(X::AbstractSingleton, Y::($ST), witness::Bool=false) = _issubset_singleton(X, Y, witness)
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
    result = element(S1) == element(S2) # TODO use _isapprox or isequivalent
    if witness
        N = promote_type(eltype(S1), eltype(S2))
        return (result, result ? N[] : element(S1))
    else
        return result
    end
end


# --- Ball2 ---


"""
    ⊆(B1::Ball2, B2::Ball2{N}, [witness]::Bool=false
     ) where {N<:AbstractFloat}

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
    if witness
        if result
            N = promote_type(eltype(B1), eltype(B2))
            return (result, N[])
        end
    else
        return result
    end

    # compute a witness 'v'
    v = B1.center .+ B1.radius * (B1.center .- B2.center)
    return (false, v)
end


# --- Ball2/Ballp ---


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
    result = B.center == element(S) && B.radius == 0 # TODO use _isapprox
    if witness
        if result
            N = promote_type(eltype(B), eltype(S))
            return (result, N[])
        end
    else
        return result
    end

    # compute a witness 'p' in the difference
    if B.center != element(S)
        p = B.center
    else
        p = copy(B.center)
        p[1] += B.radius
    end
    return (false, p)
end


# --- LineSegment ---


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
    return _issubset_line_segment(L, S)
end

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

# fix ambiguity
⊆(L::LineSegment, H::AbstractHyperrectangle, witness::Bool=false) = _issubset_line_segment(L, H, witness)

# --- Interval ---


"""
    ⊆(x::Interval, y::Interval, [witness]::Bool=false)

Check whether an interval is contained in another interval.

### Input

- `x`       -- interval
- `y`       -- interval
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

`true` iff ``x ⊆ y``.
"""
function ⊆(x::Interval, y::Interval, witness::Bool=false)
    witness && raise(ValueError("witness production is not supported yet"))
    return x.dat ⊆ y.dat
end


# --- EmptySet ---


"""
    ⊆(∅::EmptySet, X::LazySet, witness::Bool=false)

Check whether an empty set is contained in another set.

### Input

- `∅`       -- empty set
- `X`       -- another set
- `witness` -- (optional, default: `false`) compute a witness if activated
               (ignored, just kept for interface reasons)

### Output

`true`.
"""
function ⊆(∅::EmptySet, X::LazySet, witness::Bool=false)
    return _issubset_emptyset(∅, X, witness)
end

function _issubset_emptyset(∅, X, witness)
    N = promote_type(eltype(∅), eltype(X))
    return witness ? (true, N[]) : true
end

# disambiguations
⊆(∅::EmptySet, P::AbstractPolyhedron, witness::Bool=false) = _issubset_emptyset(∅, P, witness)
⊆(∅::EmptySet, H::AbstractHyperrectangle, witness::Bool=false) = _issubset_emptyset(∅, H, witness)

"""
    ⊆(X::LazySet, ∅::EmptySet, [witness]::Bool=false)

Check whether a set is contained in an empty set.

### Input

- `X`       -- another set
- `∅`       -- empty set
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

function _issubset_in_emptyset(X, ∅, witness)
    if isempty(X)
        N = promote_type(eltype(∅), eltype(X))
        return witness ? (true, N[]) : true
    else
        return witness ? (false, an_element(X)) : false
    end
end

# disambiguations
⊆(X::AbstractPolytope, ∅::EmptySet, witness::Bool=false) = _issubset_in_emptyset(X, ∅, witness)

function ⊆(X::AbstractSingleton, ::EmptySet, witness::Bool=false)
    return witness ? (false, an_element(X)) : false
end

function ⊆(X::LineSegment, ::EmptySet, witness::Bool=false)
    return witness ? (false, an_element(X)) : false
end

function ⊆(∅₁::EmptySet, ∅₂::EmptySet, witness::Bool=false)
    N = promote_type(eltype(∅₁), eltype(∅₂))
    return witness ? (true, N[]) : true
end


# --- UnionSet ---


"""
    ⊆(cup::UnionSet, X::LazySet, [witness]::Bool=false)

Check whether a union of two convex sets is contained in another set.

### Input

- `cup`     -- union of two convex sets
- `X`       -- another set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``\\text{cup} ⊆ X``
* If `witness` option is activated:
  * `(true, [])` iff ``\\text{cup} ⊆ X``
  * `(false, v)` iff ``\\text{cup} \\not\\subseteq X`` and
    ``v ∈ \\text{cup} \\setminus X``
"""
function ⊆(cup::UnionSet, X::LazySet, witness::Bool=false)
    return ⊆(UnionSetArray([cup.X, cup.Y]), X, witness) # TODO implement here
end

"""
    ⊆(cup::UnionSetArray, X::LazySet, [witness]::Bool=false)

Check whether a union of a finite number of convex sets is contained in another
set.

### Input

- `cup`     -- union of a finite number of convex sets
- `X`       -- another set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``\\text{cup} ⊆ X``
* If `witness` option is activated:
  * `(true, [])` iff ``\\text{cup} ⊆ X``
  * `(false, v)` iff ``\\text{cup} \\not\\subseteq X`` and
    ``v ∈ \\text{cup} \\setminus X``
"""
function ⊆(cup::UnionSetArray, X::LazySet, witness::Bool=false)
    result = true
    N = promote_type(eltype(cup), eltype(X))
    w = N[]
    for Y in array(cup)
        if witness
            result, w = ⊆(Y, X, witness)
        else
            result = ⊆(Y, X, witness)
        end
        if !result
            break
        end
    end
    return witness ? (result, w) : result
end


# --- Universe ---


"""
    ⊆(X::LazySet, U::Universe, [witness]::Bool=false)

Check whether a convex set is contained in a universe.

### Input

- `U`       -- universe
- `X`       -- convex set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true`
* If `witness` option is activated: `(true, [])`
"""
function ⊆(X::LazySet, U::Universe, witness::Bool=false)
    return _issubset_universe(X, U, witness)
end

function _issubset_universe(X, U, witness)
    N = promote_type(eltype(X), eltype(U))
    witness ? (true, N[]) : true
end

# disambiguations
for ST in [AbstractPolytope, AbstractHyperrectangle, AbstractSingleton, LineSegment, EmptySet]
    @eval ⊆(X::($ST), U::Universe, witness::Bool=false) = _issubset_universe(X, U, witness)
end

"""
    ⊆(U::Universe, X::LazySet, [witness]::Bool=false)

Check whether a universe is contained in another convex set, and otherwise
optionally compute a witness.

### Input

- `U`       -- universe
- `X`       -- convex set
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
function ⊆(U::Universe, X::LazySet, witness::Bool=false)
    return isuniversal(X, witness)
end

# disambiguations
function ⊆(U1::Universe, U2::Universe, witness::Bool=false)
    N = promote_type(eltype(U1), eltype(U2))
    return witness ? (true, N[]) : true
end

# disambiguations
for ST in [AbstractPolyhedron, AbstractPolytope, AbstractHyperrectangle, AbstractSingleton, EmptySet]
    @eval ⊆(U::Universe, X::($ST), witness::Bool=false) = isuniversal(X, witness)
end

# --- Complement ---


"""
    ⊆(X::LazySet, C::Complement, [witness]::Bool=false)

Check whether a convex set is contained in the complement of another convex set,
and otherwise optionally compute a witness.

### Input

- `X`       -- convex set
- `C`       -- complement of a convex set
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


# --- CartesianProduct ---


"""
    ⊆(X::CartesianProduct, Y::CartesianProduct, [witness]::Bool=false;
      check_block_equality::Bool=true)

Check whether a Cartesian product of two convex sets is contained in another
Cartesian product of two convex sets, and otherwise optionally compute a
witness.

### Input

- `X`       -- Cartesian product of two convex sets
- `Y`       -- Cartesian product of two convex sets
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

We check for inclusion for each block of the Cartesian products.

For witness production, we obtain a witness in one of the blocks.
We then construct a high-dimensional witness by obtaining any point in the other
blocks (using `an_element`) and concatenating these points.
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

    N = promote_type(eltype(X), eltype(Y))
    return witness ? (true, N[]) : true
end

"""
    ⊆(X::CartesianProductArray, Y::CartesianProductArray, [witness]::Bool=false;
      check_block_equality::Bool=true)

Check whether a Cartesian product of finitely many convex sets is contained in
another Cartesian product of finitely many convex sets, and otherwise optionally
compute a witness.

### Input

- `X`       -- Cartesian product of finitely many convex sets
- `Y`       -- Cartesian product of finitely many convex sets
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

We check for inclusion for each block of the Cartesian products.

For witness production, we obtain a witness in one of the blocks.
We then construct a high-dimensional witness by obtaining any point in the other
blocks (using `an_element`) and concatenating these points.
"""
function ⊆(X::CartesianProductArray, Y::CartesianProductArray, witness::Bool=false;
           check_block_equality::Bool=true)
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


# --- AbstractZonotope ---


"""
    ⊆(Z::AbstractZonotope, H::AbstractHyperrectangle, [witness]::Bool=false)

Check whether a zonotopic set is contained in a hyperrectangular set.

### Input

- `Z`       -- inner zonotopic set
- `H`       -- outer hyperrectangular set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

`true` iff ``Z ⊆ H`` otherwise `false`

### Algorithm

Algorithm based on Lemma 3.1 of [1]

[1] Mitchell, I. M., Budzis, J., & Bolyachevets, A. (2019, April). Invariant,
 viability and discriminating kernel under-approximation via zonotope scaling.
 In Proceedings of the 22nd ACM International Conference on Hybrid Systems:
 Computation and Control (pp. 268-269).
"""
function ⊆(Z::AbstractZonotope, H::AbstractHyperrectangle, witness::Bool=false)
    witness && raise(ValueError("witness production is not supported yet"))

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
