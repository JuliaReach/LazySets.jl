import Base.issubset

export is_subset


# --- AbstractHyperrectangle ---


"""
    is_subset(S::LazySet{N}, H::AbstractHyperrectangle{N}, [witness]::Bool=false
             )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a convex set is contained in a hyperrectangle, and if not,
optionally compute a witness.

### Input

- `S` -- inner convex set
- `H` -- outer hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``S ⊆ H``
  * `(false, v)` iff ``S \\not\\subseteq H`` and ``v ∈ S \\setminus H``

### Algorithm

``S ⊆ H`` iff ``\\operatorname{ihull}(S) ⊆ H``, where  ``\\operatorname{ihull}``
is the interval hull operator.
"""
function is_subset(S::LazySet{N},
                   H::AbstractHyperrectangle{N},
                   witness::Bool=false
                  )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return is_subset(Approximations.interval_hull(S), H, witness)
end


"""
    is_subset(P::AbstractPolytope{N},
              H::AbstractHyperrectangle,
              [witness]::Bool=false
             )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a polytope is contained in a hyperrectangle, and if not,
optionally compute a witness.

### Input

- `P` -- inner polytope
- `H` -- outer hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``P ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``P ⊆ H``
  * `(false, v)` iff ``P \\not\\subseteq H`` and ``v ∈ P \\setminus H``

### Notes

This copy-pasted method just exists to avoid method ambiguities.

### Algorithm

Since ``H`` is convex, ``P ⊆ H`` iff ``v_i ∈ H`` for all vertices ``v_i`` of
``P``.
"""
function is_subset(P::AbstractPolytope{N},
                   H::AbstractHyperrectangle,
                   witness::Bool=false
                  )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    @assert dim(P) == dim(H)

    for v in vertices_list(P)
        if !∈(v, H)
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
    is_subset(H1::AbstractHyperrectangle{N},
              H2::AbstractHyperrectangle{N},
              [witness]::Bool=false
             )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a given hyperrectangle is contained in another hyperrectangle, and
if not, optionally compute a witness.

### Input

- `H1` -- inner hyperrectangle
- `H2` -- outer hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``H1 ⊆ H2``
* If `witness` option is activated:
  * `(true, [])` iff ``H1 ⊆ H2``
  * `(false, v)` iff ``H1 \\not\\subseteq H2`` and ``v ∈ H1 \\setminus H2``

### Algorithm

``H1 ⊆ H2`` iff ``c_1 + r_1 ≤ c_2 + r_2 ∧ c_1 - r_1 ≥ c_2 - r_2`` iff
``r_1 - r_2 ≤ c_1 - c_2 ≤ -(r_1 - r_2)``, where ``≤`` is taken component-wise.
"""
function is_subset(H1::AbstractHyperrectangle{N},
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
    is_subset(P::AbstractPolytope{N}, S::LazySet{N}, [witness]::Bool=false
             )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a polytope is contained in a convex set, and if not, optionally
compute a witness.

### Input

- `P` -- inner polytope
- `S` -- outer convex set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``P ⊆ S``
* If `witness` option is activated:
  * `(true, [])` iff ``P ⊆ S``
  * `(false, v)` iff ``P \\not\\subseteq S`` and ``v ∈ P \\setminus S``

### Algorithm

Since ``S`` is convex, ``P ⊆ S`` iff ``v_i ∈ S`` for all vertices ``v_i`` of
``P``.
"""
function is_subset(P::AbstractPolytope{N}, S::LazySet{N}, witness::Bool=false
                  )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    @assert dim(P) == dim(S)

    for v in vertices_list(P)
        if !∈(v, S)
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


# --- AbstractSingleton ---


"""
    is_subset(S::AbstractSingleton{N}, set::LazySet{N}, [witness]::Bool=false
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
  * `(false, v)` iff ``S \\not\\subseteq \\text{set}`` and
    ``v ∈ S \\setminus \\text{set}``
"""
function is_subset(S::AbstractSingleton{N}, set::LazySet{N}, witness::Bool=false
                  )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    result = ∈(element(S), set)
    if witness
        return (result, result ? N[] : element(S))
    else
        return result
    end
end


"""
    is_subset(S::AbstractSingleton{N},
              H::AbstractHyperrectangle{N},
              [witness]::Bool=false
             )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a given set with a single value is contained in a hyperrectangle,
and if not, optionally compute a witness.

### Input

- `S` -- inner set with a single value
- `H` -- outer hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ⊆ H``
* If `witness` option is activated:
  * `(true, [])` iff ``S ⊆ H``
  * `(false, v)` iff ``S \\not\\subseteq H`` and ``v ∈ S \\setminus H``

### Notes

This copy-pasted method just exists to avoid method ambiguities.
"""
function is_subset(S::AbstractSingleton{N},
                   H::AbstractHyperrectangle{N},
                   witness::Bool=false
                  )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    result = ∈(element(S), H)
    if witness
        return (result, result ? N[] : element(S))
    else
        return result
    end
end


"""
    is_subset(S1::AbstractSingleton{N},
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
  * `(false, v)` iff ``S1 \\not\\subseteq S2`` and ``v ∈ S1 \\setminus S2``
"""
function is_subset(S1::AbstractSingleton{N},
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
    is_subset(B1::Ball2{N}, B2::Ball2{N}, [witness]::Bool=false
             )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
  * `(false, v)` iff ``B1 \\not\\subseteq B2`` and ``v ∈ B1 \\setminus B2``

### Algorithm

``B1 ⊆ B2`` iff ``‖ c_1 - c_2 ‖_2 + r_1 ≤ r_2``
"""
function is_subset(B1::Ball2{N}, B2::Ball2{N}, witness::Bool=false
                  )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
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
    is_subset(B::Union{Ball2{N}, Ballp{N}},
              S::AbstractSingleton{N},
              [witness]::Bool=false
             )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

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
  * `(false, v)` iff ``B \\not\\subseteq S`` and ``v ∈ B \\setminus S``
"""
function is_subset(B::Union{Ball2{N}, Ballp{N}},
                   S::AbstractSingleton{N},
                   witness::Bool=false
                  )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
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

# TODO commented until AbstractConvexSet is implemented
# """
#     is_subset(L::LineSegment{N},
#               S::AbstractConvexSet{N},
#               witness::Bool=false
#              )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
#
# Check whether a line segment is contained in a convex set, and if not,
# optionally compute a witness.
#
# ### Input
#
# - `L` -- inner line segment
# - `S` -- outer convex set
# - `witness` -- (optional, default: `false`) compute a witness if activated
#
# ### Output
#
# * If `witness` option is deactivated: `true` iff ``L ⊆ S``
# * If `witness` option is activated:
#   * `(true, [])` iff ``L ⊆ S``
#   * `(false, v)` iff ``L \\not\\subseteq S`` and ``v ∈ L \\setminus S``
#
# ### Algorithm
#
# Since ``S`` is convex, ``L ⊆ S`` iff ``p ∈ S`` and ``q ∈ S``, where ``p, q`` are
# the end points of ``L``.
# """
# function is_subset(L::LineSegment{N},
#                    S::AbstractConvexSet{N},
#                    witness::Bool=false
#                   )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
#     p_in_S = ∈(L.p, S)
#     result = p_in_S && ∈(L.q, S)
#     if !witness
#         return result
#     elseif result
#         return (result, N[])
#     else
#         return (result, p_in_S ? L.q : L.p)
#     end
# end


# --- alias ---


"""
    ⊆

Alias for `is_subset`.

### Notes

`⊆` is a Julia-internal function which is defined for every type combination,
but crashes with a cryptic error message if it is not implemented:

    `MethodError: no method matching start(::FIRST_SET_TYPE)`
"""
@inline function issubset(S1::LazySet{N}, S2::LazySet{N}, witness::Bool=false
                         )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return is_subset(S1, S2, witness)
end
