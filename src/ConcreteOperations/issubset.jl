import Base.issubset

# this operation is forbidden, but it is a common error so we give a detailed error message
function ⊆(::AbstractVector{N}, ::LazySet{N})::Bool where {N<:Real}
    error("cannot make an inclusion check if the left-hand side " *
          "is a vector; either wrap it as a set with one element, as in " *
          "`Singleton(v) ⊆ X`, or check for set membership, as in `v ∈ X` " *
          "(they behave equivalently although the implementations may differ)")
end


# --- fall-back with constraints_list of rhs ---

"""
    ⊆(X::LazySet{N}, P::LazySet{N}, [witness]::Bool=false) where {N<:Real}

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
function ⊆(X::LazySet{N}, P::LazySet{N}, witness::Bool=false) where {N<:Real}
    if applicable(constraints_list, P)
        return _issubset_constraints_list(X, P, witness)
    else
        error("an inclusion check for the given combination of set types is " *
              "not available")
    end
end


# --- AbstractHyperrectangle ---


"""
    ⊆(S::LazySet{N}, H::AbstractHyperrectangle{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function ⊆(S::LazySet{N}, H::AbstractHyperrectangle{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return ⊆(Approximations.interval_hull(S), H, witness)
end

"""
    ⊆(P::AbstractPolytope{N}, H::AbstractHyperrectangle, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a polytope is contained in a hyperrectangular set, and if not,
optionally compute a witness.

### Input

- `P` -- inner polytope
- `H` -- outer hyperrectangular set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``P ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``P ⊆ H``
  * `(false, v)` iff ``P ⊈ H`` and ``v ∈ P \\setminus H``

### Notes

This copy-pasted method just exists to avoid method ambiguities.

### Algorithm

Since ``H`` is convex, ``P ⊆ H`` iff ``v_i ∈ H`` for all vertices ``v_i`` of
``P``.
"""
function ⊆(P::AbstractPolytope{N},
           H::AbstractHyperrectangle,
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    @assert dim(P) == dim(H)

    return _issubset_vertices_list(P, H, witness)
end


"""
    ⊆(H1::AbstractHyperrectangle{N},
      H2::AbstractHyperrectangle{N},
      [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function ⊆(H1::AbstractHyperrectangle{N},
           H2::AbstractHyperrectangle{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    @assert dim(H1) == dim(H2)

    for i in 1:dim(H1)
        c_dist = center(H1)[i] - center(H2)[i]
        r_dist = radius_hyperrectangle(H1, i) - radius_hyperrectangle(H2, i)
        if -r_dist < c_dist || c_dist < r_dist
            if witness
                # compute a witness 'p' in the difference
                p = copy(center(H1))
                if c_dist >= 0
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
    ⊆(P::AbstractPolytope{N}, S::LazySet{N}, [witness]::Bool=false;
      algorithm="constraints")::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a polytope is contained in a convex set, and if not, optionally
compute a witness.

### Input

- `P` -- inner polytope
- `S` -- outer convex set
- `witness`   -- (optional, default: `false`) compute a witness if activated
- `algorithm` -- (optional, default: `"constraints"`) algorithm for the inclusion
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
function ⊆(P::AbstractPolytope{N}, S::LazySet{N}, witness::Bool=false;
           algorithm="constraints")::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    @assert dim(P) == dim(S)

    if algorithm == "constraints"
        return _issubset_constraints_list(P, S, witness)
    elseif algorithm == "vertices"
        return _issubset_vertices_list(P, S, witness)
    else
        error("algorithm $algorithm unknown")
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
    ⊆(X::LazySet{N},
      P::AbstractPolyhedron{N},
      witness::Bool=false
     ) where {N<:Real}

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
function ⊆(X::LazySet{N}, P::AbstractPolyhedron{N}, witness::Bool=false
          ) where {N<:Real}
    return _issubset_constraints_list(X, P, witness)
end

# for documentation see
# ⊆(X::LazySet{N}, P::AbstractPolyhedron{N}, witness::Bool=false) where {N<:Real}
function _issubset_constraints_list(S::LazySet{N}, P::LazySet{N},
                                    witness::Bool=false) where {N<:Real}
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
    return witness ? (true, N[]) : true
end

# disambiguation
function ⊆(P1::AbstractPolytope{N},
           P2::AbstractPolyhedron{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return invoke(⊆, Tuple{LazySet{N}, typeof(P2), Bool}, P1, P2, witness)
end
function ⊆(H::AbstractHyperrectangle{N},
           P::AbstractPolyhedron{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return invoke(⊆, Tuple{LazySet{N}, typeof(P), Bool}, H, P, witness)
end
function ⊆(S::AbstractSingleton{N},
           P::AbstractPolyhedron{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return invoke(⊆, Tuple{LazySet{N}, typeof(P), Bool}, S, P, witness)
end
function ⊆(L::LineSegment{N},
           P::AbstractPolyhedron{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return invoke(⊆, Tuple{LazySet{N}, typeof(P), Bool}, L, P, witness)
end
function ⊆(P::AbstractPolytope{N},
           H::AbstractHyperrectangle{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return invoke(⊆,
                  Tuple{LazySet{N}, AbstractHyperrectangle{N}, Bool},
                  P, H, witness)
end


# --- AbstractSingleton ---


"""
    ⊆(S::AbstractSingleton{N}, set::LazySet{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a given set with a single value is contained in a convex set, and
if not, optionally compute a witness.

### Input

- `S`   -- inner set with a single value
- `set` -- outer convex set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ⊆ \\text{set}``
* If `witness` option is activated:
  * `(true, [])` iff ``S ⊆ \\text{set}``
  * `(false, v)` iff ``S ⊈ \\text{set}`` and
    ``v ∈ S \\setminus \\text{set}``
"""
function ⊆(S::AbstractSingleton{N}, set::LazySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    result = element(S) ∈ set
    if witness
        return (result, result ? N[] : element(S))
    else
        return result
    end
end


"""
    ⊆(S::AbstractSingleton{N},
      H::AbstractHyperrectangle{N},
      [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a given set with a single value is contained in a hyperrectangular
set, and if not, optionally compute a witness.

### Input

- `S` -- inner set with a single value
- `H` -- outer hyperrectangular set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``S ⊆ H``
  * `(false, v)` iff ``S ⊈ H`` and ``v ∈ S \\setminus H``

### Notes

This copy-pasted method just exists to avoid method ambiguities.
"""
function ⊆(S::AbstractSingleton{N},
           H::AbstractHyperrectangle{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    result = element(S) ∈ H
    if witness
        return (result, result ? N[] : element(S))
    else
        return result
    end
end


"""
    ⊆(S1::AbstractSingleton{N},
      S2::AbstractSingleton{N},
      [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function ⊆(S1::AbstractSingleton{N},
           S2::AbstractSingleton{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    result = element(S1) == element(S2)
    if witness
        return (result, result ? N[] : element(S1))
    else
        return result
    end
end


# --- Ball2 ---


"""
    ⊆(B1::Ball2{N}, B2::Ball2{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:AbstractFloat}

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
function ⊆(B1::Ball2{N}, B2::Ball2{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:AbstractFloat}
    result = norm(B1.center - B2.center, 2) + B1.radius <= B2.radius
    if witness
        if result
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
    ⊆(B::Union{Ball2{N}, Ballp{N}},
      S::AbstractSingleton{N},
      [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:AbstractFloat}

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
function ⊆(B::Union{Ball2{N}, Ballp{N}},
           S::AbstractSingleton{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:AbstractFloat}
    result = B.center == element(S) && B.radius == 0
    if witness
        if result
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
    ⊆(L::LineSegment{N}, S::LazySet{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function ⊆(L::LineSegment{N}, S::LazySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    p_in_S = L.p ∈ S
    result = p_in_S && L.q ∈ S
    if !witness
        return result
    elseif result
        return (result, N[])
    else
        return (result, p_in_S ? L.q : L.p)
    end
end


"""
    ⊆(L::LineSegment{N}, H::AbstractHyperrectangle{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a line segment is contained in a hyperrectangular set, and if not,
optionally compute a witness.

### Input

- `L` -- inner line segment
- `H` -- outer hyperrectangular set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``L ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``L ⊆ H``
  * `(false, v)` iff ``L ⊈ H`` and ``v ∈ L \\setminus H``

### Notes

This copy-pasted method just exists to avoid method ambiguities.

### Algorithm

Since ``H`` is convex, ``L ⊆ H`` iff ``p ∈ H`` and ``q ∈ H``, where ``p, q`` are
the end points of ``L``.
"""
function ⊆(L::LineSegment{N}, H::AbstractHyperrectangle{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    p_in_H = L.p ∈ H
    result = p_in_H && L.q ∈ H
    if !witness
        return result
    elseif result
        return (result, N[])
    else
        return (result, p_in_H ? L.q : L.p)
    end
end


# --- Interval ---


"""
    ⊆(x::Interval, y::Interval)

Check whether an interval is contained in another interval.

### Input

- `x` -- interval
- `y` -- interval

### Output

`true` iff ``x ⊆ y``.
"""
function ⊆(x::Interval, y::Interval)
    return x.dat ⊆ y.dat
end


# --- EmptySet ---


"""
    ⊆(∅::EmptySet{N}, X::LazySet{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether an empty set is contained in another set.

### Input

- `∅`       -- empty set
- `X`       -- another set
- `witness` -- (optional, default: `false`) compute a witness if activated
               (ignored, just kept for interface reasons)

### Output

`true`.
"""
function ⊆(∅::EmptySet{N}, X::LazySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end

# disambiguation
function ⊆(∅::EmptySet{N},
           ::AbstractPolyhedron{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end
function ⊆(∅::EmptySet{N},
           ::AbstractHyperrectangle{N},
           witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end

"""
    ⊆(X::LazySet{N}, ∅::EmptySet{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function ⊆(X::LazySet{N}, ∅::EmptySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    if isempty(X)
        return witness ? (true, N[]) : true
    else
        return witness ? (false, an_element(X)) : false
    end
end

# disambiguation
function ⊆(X::AbstractPolytope{N}, ::EmptySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    if isempty(X)
        return witness ? (true, N[]) : true
    else
        return witness ? (false, an_element(X)) : false
    end
end
function ⊆(X::AbstractSingleton{N}, ::EmptySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (false, an_element(X)) : false
end
function ⊆(X::LineSegment{N}, ::EmptySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (false, an_element(X)) : false
end
function ⊆(::EmptySet{N}, ::EmptySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end


# --- UnionSet ---


"""
    ⊆(cup::UnionSet{N}, X::LazySet{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function ⊆(cup::UnionSet{N}, X::LazySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return ⊆(UnionSetArray([cup.X, cup.Y]), X, witness)
end

"""
    ⊆(cup::UnionSetArray{N}, X::LazySet{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function ⊆(cup::UnionSetArray{N}, X::LazySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    result = true
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
    ⊆(X::LazySet{N}, U::Universe{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a convex set is contained in a universe.

### Input

- `U`       -- universe
- `X`       -- convex set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true`
* If `witness` option is activated: `(true, [])`
"""
function ⊆(X::LazySet{N}, U::Universe{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end

# disambiguation
function ⊆(P::AbstractPolytope{N}, ::Universe{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end
function ⊆(P::AbstractHyperrectangle{N}, ::Universe{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end
function ⊆(S::AbstractSingleton{N}, ::Universe{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end
function ⊆(::LineSegment{N}, ::Universe{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end
function ⊆(::EmptySet{N}, ::Universe{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end

"""
    ⊆(U::Universe{N}, X::LazySet{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function ⊆(U::Universe{N}, X::LazySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return isuniversal(X, witness)
end

# disambiguation
function ⊆(::Universe{N}, ::Universe{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return witness ? (true, N[]) : true
end
function ⊆(::Universe{N}, P::AbstractPolyhedron{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return isuniversal(P, witness)
end
function ⊆(::Universe{N}, P::AbstractPolytope{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return isuniversal(P, witness)
end
function ⊆(::Universe{N}, H::AbstractHyperrectangle{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return isuniversal(H, witness)
end
function ⊆(::Universe{N}, S::AbstractSingleton{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return isuniversal(S, witness)
end
function ⊆(::Universe{N}, ::EmptySet{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return isuniversal(P, witness)
end


# --- Complement ---


"""
    ⊆(X::LazySet{N}, C::Complement{N}, [witness]::Bool=false
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
function ⊆(X::LazySet{N}, C::Complement{N}, witness::Bool=false
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return isdisjoint(X, C.X, witness)
end


# --- CartesianProduct ---


"""
    ⊆(X::CartesianProductArray{N}, Y::CartesianProductArray{N},
      witness::Bool=false; check_block_equality::Bool=true
     )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
Depending on the value of `check_block_equality`, we check this property.

### Algorithm

We check for inclusion for each block of the Cartesian products.

For witness production, we obtain a witness in one of the blocks.
We then construct a high-dimensional witness by obtaining any point in the other
blocks (using `an_element`) and concatenating these points.
"""
function ⊆(X::CartesianProductArray{N},
           Y::CartesianProductArray{N},
           witness::Bool=false;
           check_block_equality::Bool=true
          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    aX = array(X)
    aY = array(Y)
    if check_block_equality && !same_block_structure(aX, aY)
        throw(ArgumentError("this inclusion check requires Cartesian products" *
                            "with the same block structure"))
    end

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
