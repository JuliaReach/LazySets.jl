export is_intersection_empty,
       isdisjoint


# --- AbstractHyperrectangle ---


"""
    is_intersection_empty(H1::AbstractHyperrectangle{N},
                          H2::AbstractHyperrectangle{N},
                          witness::Bool=false
                         )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether two hyperrectangles do not intersect, and otherwise optionally
compute a witness.

### Input

- `H1` -- first hyperrectangle
- `H2` -- second hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``H1 ∩ H2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``H1 ∩ H2 = ∅``
  * `(false, v)` iff ``H1 ∩ H2 ≠ ∅`` and ``v ∈ H1 ∩ H2``

### Algorithm

``H1 ∩ H2 ≠ ∅`` iff ``|c_2 - c_1| ≤ r_1 + r_2``, where ``≤`` is taken
component-wise.

A witness is computed by starting in one center and moving toward the other
center for as long as the minimum of the radius and the center distance.
In other words, the witness is the point in `H1` that is closest to the center
of `H2`.
"""
function is_intersection_empty(H1::AbstractHyperrectangle{N},
                               H2::AbstractHyperrectangle{N},
                               witness::Bool=false
                              )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    empty_intersection = false
    center_diff = center(H2) - center(H1)
    for i in eachindex(center_diff)
        if abs(center_diff[i]) >
                radius_hyperrectangle(H1, i) + radius_hyperrectangle(H2, i)
            empty_intersection = true
            break
        end
    end

    if !witness
        return empty_intersection
    elseif empty_intersection
        return (true, N[])
    end

    # compute a witness 'v' in the intersection
    v = copy(center(H1))
    c2 = center(H2)
    for i in eachindex(center_diff)
        if v[i] <= c2[i]
            # second center is right of first center
            v[i] += min(radius_hyperrectangle(H1, i), center_diff[i])
        else
            # second center is left of first center
            v[i] -= max(radius_hyperrectangle(H1, i), -center_diff[i])
        end
    end
    return (false, v)
end


# --- AbstractSingleton ---


"""
    is_intersection_empty(S::AbstractSingleton{N},
                          set::LazySet{N},
                          witness::Bool=false
                         )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a singleton and a convex set do not intersect, and otherwise
optionally compute a witness.

### Input

- `S`   -- singleton
- `set` -- convex set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ∩ \\operatorname{set} = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``S ∩ \\operatorname{set} = ∅``
  * `(false, v)` iff ``S ∩ \\operatorname{set} ≠ ∅`` and
    `v` = `element(S)` ``∈ S ∩ \\operatorname{set}``

### Algorithm

``S ∩ \\operatorname{set} = ∅`` iff `element(S)` ``∉ \\operatorname{set}``.
"""
function is_intersection_empty(S::AbstractSingleton{N},
                               set::LazySet{N},
                               witness::Bool=false
                              )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    empty_intersection = element(S) ∉ set
    if witness
        return (empty_intersection, empty_intersection ? N[] : element(S))
    else
        return empty_intersection
    end
end


"""
    is_intersection_empty(set::LazySet{N},
                          S::AbstractSingleton{N},
                          witness::Bool=false
                         )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a convex set and a singleton do not intersect, and otherwise
optionally compute a witness.

### Input

- `set` -- convex set
- `S`   -- singleton
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ∩ \\operatorname{set} = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``S ∩ \\operatorname{set} = ∅``
  * `(false, v)` iff ``S ∩ \\operatorname{set} ≠ ∅`` and
    `v` = `element(S)` ``∈ S ∩ \\operatorname{set}``

### Algorithm

``S ∩ \\operatorname{set} = ∅`` iff `element(S)` ``∉ \\operatorname{set}``.
"""
function is_intersection_empty(set::LazySet{N},
                               S::AbstractSingleton{N},
                               witness::Bool=false
                              )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    return is_intersection_empty(S, set, witness)
end


"""
    is_intersection_empty(S1::AbstractSingleton{N},
                          S2::AbstractSingleton{N},
                          witness::Bool=false
                         )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether two singletons do not intersect, and otherwise optionally compute
a witness.

### Input

- `S1` -- first singleton
- `S2` -- second singleton
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S1 ∩ S2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``S1 ∩ S2 = ∅``
  * `(false, v)` iff ``S1 ∩ S2 ≠ ∅`` and `v` = `element(S1)` ``∈ S1 ∩ S2``

### Algorithm

``S1 ∩ S2 = ∅`` iff ``S1 ≠ S2``.
"""
function is_intersection_empty(S1::AbstractSingleton{N},
                               S2::AbstractSingleton{N},
                               witness::Bool=false
                              )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    empty_intersection = element(S1) != element(S2)
    if witness
        return (empty_intersection, empty_intersection ? N[] : element(S1))
    else
        return empty_intersection
    end
end


"""
    is_intersection_empty(H::AbstractHyperrectangle{N},
                          S::AbstractSingleton{N},
                          witness::Bool=false
                         )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a hyperrectangle and a singleton do not intersect, and otherwise
optionally compute a witness.

### Input

- `H` -- hyperrectangle
- `S` -- singleton
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``H ∩ S = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``H ∩ S = ∅``
  * `(false, v)` iff ``H ∩ S ≠ ∅`` and `v` = `element(S)` ``∈ H ∩ S``

### Algorithm

``H ∩ S = ∅`` iff `element(S)` ``∉ H``.
"""
function is_intersection_empty(H::AbstractHyperrectangle{N},
                               S::AbstractSingleton{N},
                               witness::Bool=false
                              )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    empty_intersection = element(S) ∉ H
    if witness
        return (empty_intersection, empty_intersection ? N[] : element(S))
    else
        return empty_intersection
    end
end


"""
    is_intersection_empty(S::AbstractSingleton{N},
                          H::AbstractHyperrectangle{N},
                          witness::Bool=false
                         )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether a singleton and a hyperrectangle do not intersect, and otherwise
optionally compute a witness.

### Input

- `S` -- singleton
- `H` -- hyperrectangle
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``H ∩ S = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``H ∩ S = ∅``
  * `(false, v)` iff ``H ∩ S ≠ ∅`` and `v` = `element(S)` ``∈ H ∩ S``

### Algorithm

``S ∩ H = ∅`` iff `element(S)` ``∉ H``.
"""
function is_intersection_empty(S::AbstractSingleton{N},
                               H::AbstractHyperrectangle{N},
                               witness::Bool=false
                              )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    return is_intersection_empty(H, S, witness)
end


# --- Ball2 ---


"""
    is_intersection_empty(B1::Ball2{N},
                          B2::Ball2{N},
                          witness::Bool=false
                         )::Union{Bool, Tuple{Bool,Vector{N}}} where {N<:Real}

Check whether two balls in the 2-norm do not intersect, and otherwise optionally
compute a witness.

### Input

- `B1` -- first ball in the 2-norm
- `B2` -- second ball in the 2-norm
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``B1 ∩ B2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``B1 ∩ B2 = ∅``
  * `(false, v)` iff ``B1 ∩ B2 ≠ ∅`` and ``v ∈ B1 ∩ B2``

### Algorithm

``B1 ∩ B2 = ∅`` iff ``‖ c_2 - c_1 ‖_2 > r_1 + r_2``.

A witness is computed depending on the smaller/bigger ball (to break ties,
choose `B1` for the smaller ball) as follows.
- If the smaller ball's center is contained in the bigger ball, we return it.
- Otherwise start in the smaller ball's center and move toward the other center
  until hitting the smaller ball's border.
  In other words, the witness is the point in the smaller ball that is closest
  to the center of the bigger ball.
"""
function is_intersection_empty(B1::Ball2{N},
                               B2::Ball2{N},
                               witness::Bool=false
                              )::Union{Bool, Tuple{Bool,Vector{N}}} where {N<:Real}
    center_diff_normed = norm(center(B2) - center(B1), 2)
    empty_intersection = center_diff_normed > B1.radius + B2.radius

    if !witness
        return empty_intersection
    elseif empty_intersection
        return (true, N[])
    end

    # compute a witness 'v' in the intersection
    if B1.radius <= B2.radius
        smaller = B1
        bigger = B2
    else
        smaller = B2
        bigger = B1
    end
    if center_diff_normed <= bigger.radius
        # smaller ball's center is contained in bigger ball
        v = smaller.center
    else
        # scale center difference with smaller ball's radius
        direction = (bigger.center - smaller.center)
        v = smaller.center + direction / center_diff_normed * smaller.radius
    end
    return (false, v)
end


# --- Zonotope ---


"""
    is_intersection_empty(Z::Zonotope{N}, H::Hyperplane{N}, witness::Bool=false
                         )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a zonotope and a hyperplane do not intersect, and otherwise
optionally compute a witness.

### Input

- `Z` -- zonotope
- `H` -- hyperplane
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``Z ∩ H = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``Z ∩ H = ∅``
  * `(false, v)` iff ``Z ∩ H ≠ ∅`` and ``v ∈ Z ∩ H``

### Algorithm

``Z ∩ H = ∅`` iff ``(b - a⋅c) ∉ \\left[ ± ∑_{i=1}^p |a⋅g_i| \\right]``,
where ``a``, ``b`` are the hyperplane coefficients, ``c`` is the zonotope's
center, and ``g_i`` are the zonotope's generators.

### Notes

Witness production is currently not supported.
"""
function is_intersection_empty(Z::Zonotope{N},
                               H::Hyperplane{N},
                               witness::Bool=false
                              )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    v = H.b - dot(H.a, Z.center)
    abs_sum = sum(abs(dot(H.a, Z.generators[:, i])) for i = 1:ngens(Z))
    empty_intersection = v < -abs_sum || v > abs_sum

    if !witness
        return empty_intersection
    elseif empty_intersection
        return (true, N[])
    else
        error("witness production is not supported yet")
    end
end

"""
    is_intersection_empty(H::Hyperplane{N}, Z::Zonotope{N}, witness::Bool=false
                         )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a hyperplane and a zonotope do not intersect, and otherwise
optionally compute a witness.

### Input

- `H` -- hyperplane
- `Z` -- zonotope
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``H ∩ Z = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``H ∩ Z = ∅``
  * `(false, v)` iff ``H ∩ Z ≠ ∅`` and ``v ∈ H ∩ Z``

### Algorithm

``H ∩ Z = ∅`` iff ``(b - a⋅c) ∉ \\left[ ± ∑_{i=1}^p |a⋅g_i| \\right]``,
where ``a``, ``b`` are the hyperplane coefficients, ``c`` is the zonotope's
center, and ``g_i`` are the zonotope's generators.

### Notes

Witness production is currently not supported.
"""
function is_intersection_empty(H::Hyperplane{N},
                               Z::Zonotope{N},
                               witness::Bool=false
                              )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return is_intersection_empty(Z, H, witness)
end

"""
    is_intersection_empty(ls1::LineSegment{N},
                          ls2::LineSegment{N},
                          witness::Bool=false
                         )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether two line segments do not intersect, and otherwise optionally
compute a witness.

### Input

- `ls1` -- first line segment
- `ls2` -- second line segment
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``ls1 ∩ ls2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``ls1 ∩ ls2 = ∅``
  * `(false, v)` iff ``ls1 ∩ ls2 ≠ ∅`` and ``v ∈ ls1 ∩ ls2``

### Algorithm

The algorithm is inspired from [here](https://stackoverflow.com/a/565282), which
again is the special 2D case of a 3D algorithm by Ronald Goldman's article on the
*Intersection of two lines in three-space* in Graphics Gems, Andrew S. (ed.), 1990.

We first check if the two line segments are parallel, and if so, if they are
collinear.
In the latter case, we check containment of any of the end points in the other
line segment.
Otherwise the lines are not parallel, so we can solve an equation of the
intersection point, if it exists.
"""
function is_intersection_empty(ls1::LineSegment{N},
                               ls2::LineSegment{N},
                               witness::Bool=false
                              )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    function cross(ls1::Vector{N}, ls2::Vector{N})::N where {N<:Real}
        return ls1[1] * ls2[2] - ls1[2] * ls2[1]
    end

    r = ls1.q - ls1.p
    if iszero(r)
        # first line segment is a point
        empty_intersection = ls1.q ∉ ls2
        if witness
            return (empty_intersection, empty_intersection ? N[] : ls1.q)
        else
            return empty_intersection
        end
    end

    s = ls2.q - ls2.p
    if iszero(s)
        # second line segment is a point
        empty_intersection = ls2.q ∉ ls1
        if witness
            return (empty_intersection, empty_intersection ? N[] : ls2.q)
        else
            return empty_intersection
        end
    end

    p1p2 = ls2.p - ls1.p
    u_numerator = cross(p1p2, r)
    u_denominator = cross(r, s)

    if u_denominator == 0
        # line segments are parallel
        if u_numerator == 0
            # line segments are collinear
            if ls1.p ∈ ls2
                empty_intersection = false
                if witness
                    v = ls1.p
                end
            elseif ls1.q ∈ ls2
                empty_intersection = false
                if witness
                    v = ls1.q
                end
            else
                empty_intersection = true
            end
        else
            # line segments are parallel and not collinear
            empty_intersection = true
        end
    else
        # line segments are not parallel
        u = u_numerator / u_denominator
        if u < 0 || u > 1
            empty_intersection = true
        else
            t = cross(p1p2, s) / u_denominator
            empty_intersection = t < 0 || t > 1
            if witness
                v = ls1.p + t * r
            end
        end
    end

    if witness
        return (empty_intersection, empty_intersection ? N[] : v)
    else
        return empty_intersection
    end
end

# --- AbstractPolytope ---

"""
    is_intersection_empty(P::AbstractPolytope{N},
                          Q::AbstractPolytope{N},
                          witness::Bool=false
                          )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether two polytopes do not intersect, and otherwise optionally
compute a witness.

### Input

- `P`       -- polytope
- `Q`       -- another polytope
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``P ∩ Q = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``P ∩ Q = ∅``
  * `(false, v)` iff ``P ∩ Q ≠ ∅`` and ``v ∈ P ∩ Q``

### Algorithm

This is a fallback implementation of the `AbstractPolytope` interface that
computes the concrete intersection, `intersection`, of the given pair of
polytopes.
If a witness is required, the first vertex of the resulting intersection
polytope is returned.
"""
function is_intersection_empty(P::AbstractPolytope{N},
                               Q::AbstractPolytope{N},
                               witness::Bool=false
                              )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    X = intersection(P, Q)
    isempty_flag = isempty(X)
    if witness
        witness_vertex = isempty_flag ? N[] : vertices_list(X)[1]
        return (isempty_flag, witness_vertex)
    else
        return isempty_flag
    end
end

# method disambiguation
function is_intersection_empty(point::AbstractSingleton{N},
                               set::AbstractPolytope{N},
                               witness::Bool=false
                              )::Union{Bool,Tuple{Bool,Vector{N}}} where {N<:Real}
    empty_intersection = element(point) ∉ set
    if witness
        return (empty_intersection, empty_intersection ? N[] : element(S))
    else
        return empty_intersection
    end
end

# symmetric function
function is_intersection_empty(set::AbstractPolytope{N},
                               point::AbstractSingleton{N},
                               witness::Bool=false
                              )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    return is_intersection_empty(point, set, witness)
end

# --- Hyperplane ---

"""
    is_intersection_empty(X::LazySet{N},
                          hp::Union{Hyperplane{N}, Line{N}},
                          [witness]::Bool=false
                         )::Union{Bool, Tuple{Bool,Vector{N}}} where N<:Real

Check whether a compact set an a hyperplane do not intersect, and otherwise
optionally compute a witness.

### Input

- `X`       -- compact set
- `hp`      -- hyperplane
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ∩ hp = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``X ∩ hp = ∅``
  * `(false, v)` iff ``X ∩ hp ≠ ∅`` and ``v ∈ X ∩ hp``

### Notes

We assume that `X` is compact.
Otherwise, the support vector queries may fail.

### Algorithm

A compact convex set intersects with a hyperplane iff the support function in
the negative resp. positive direction of the hyperplane's normal vector ``a`` is
to the left resp. right of the hyperplane's constraint ``b``:

```math
-ρ(-a) ≤ b ≤ ρ(a)
```

For witness generation, we compute a line connecting the support vectors to the
left and right, and then take the intersection of the line with the hyperplane.
We follow
[this algorithm](https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection#Algebraic_form)
for the line-hyperplane intersection.
"""
function is_intersection_empty(X::LazySet{N},
                               hp::Union{Hyperplane{N}, Line{N}},
                               witness::Bool=false
                              )::Union{Bool, Tuple{Bool,Vector{N}}} where N<:Real
    normal_hp = hp.a
    sv_left = σ(-normal_hp, X)
    if -dot(sv_left, -normal_hp) <= hp.b
        sv_right = σ(normal_hp, X)
        empty_intersection = (hp.b > dot(sv_right, normal_hp))
    else
        empty_intersection = true
    end
    if witness
        if empty_intersection
            v = N[]
        else
            point_hp = an_element(hp)
            point_line = sv_left
            dir_line = sv_right - sv_left
            d = dot((point_hp - point_line), normal_hp) /
                dot(dir_line, normal_hp)
            v = d * dir_line + point_line
        end
        return (empty_intersection, v)
    else
        return empty_intersection
    end
end

# symmetric function
function is_intersection_empty(hp::Union{Hyperplane{N}, Line{N}},
                               X::LazySet{N},
                               witness::Bool=false
                              )::Union{Bool, Tuple{Bool,Vector{N}}} where N<:Real
    return is_intersection_empty(X, hp, witness)
end

# --- HalfSpace ---

"""
    is_intersection_empty(X::LazySet{N},
                          hs::HalfSpace{N},
                          [witness]::Bool=false
                         )::Union{Bool, Tuple{Bool,Vector{N}}} where N<:Real

Check whether a compact set an a half-space do not intersect, and otherwise
optionally compute a witness.

### Input

- `X`       -- compact set
- `hs`      -- half-space
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ∩ hs = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``X ∩ hs = ∅``
  * `(false, v)` iff ``X ∩ hs ≠ ∅`` and ``v ∈ X ∩ hs``

### Notes

We assume that `X` is compact.
Otherwise, the support vector queries may fail.

### Algorithm

A compact convex set intersects with a half-space iff the support vector in
the negative direction of the half-space's normal vector ``a`` is contained in
the half-space: ``σ(-a) ∈ hs``.
The support vector is thus also a witness.

Optional keyword arguments can be passed to the `ρ` function. In particular, if
`X` is a lazy intersection, the option `upper_bound=true` may be useful.
"""
function is_intersection_empty(X::LazySet{N},
                               hs::HalfSpace{N},
                               witness::Bool=false;
                               upper_bound::Bool=false,
                               kwargs...
                               )::Union{Bool, Tuple{Bool,Vector{N}}} where N<:Real
    if !witness
        ρ_rec = upper_bound ? Approximations.ρ_upper_bound : ρ
        return -ρ_rec(-hs.a, X; kwargs...) > hs.b
    end

    # for witness production, we compute the support vector instead
    svec = σ(-hs.a, X)
    empty_intersection = svec ∉ hs
    v = empty_intersection ? N[] : svec
    return (empty_intersection, v)
end

# symmetric function
function is_intersection_empty(hs::HalfSpace{N},
                               X::LazySet{N},
                               witness::Bool=false;
                               kwargs...
                              )::Union{Bool, Tuple{Bool,Vector{N}}} where N<:Real
    return is_intersection_empty(X, hs, witness; kwargs...)
end

"""
    is_intersection_empty(X::LazySet{N},
                          P::Union{HPolytope{N}, AbstractHPolygon{N}}
                         )::Bool where {N<:Real}

Check whether a compact set and a polytope do not intersect.

### Input

- `X`  -- compact set
- `P`  -- polytope or polygon in constraint-representation

### Output

`true` iff ``X ∩ P = ∅``.

### Notes

We assume that `X` is compact. Otherwise, the support vector queries may fail.
Witness production is not supported.

### Algorithm

The algorithm relies on the intersection check between the set ``X`` and each constraint
in ``P``. It costs ``m`` support vector evaluations of ``X``, where ``m`` is
the number of constraints in ``P``.

Note that this method can be used with any set ``P`` whose constraints are known.
"""
function is_intersection_empty(X::LazySet{N},
                               P::Union{HPolytope{N}, AbstractHPolygon{N}}
                              )::Bool where N<:Real

    for Hi in constraints_list(P)
        if is_intersection_empty(X, Hi)
            return true
        end
    end
    return false
end

# symmetric function
function is_intersection_empty(P::Union{HPolytope{N}, AbstractHPolygon{N}},
                               X::LazySet{N}
                              )::Bool where N<:Real
    return is_intersection_empty(X, P)
end

# --- alias ---


"""
    isdisjoint(X, Y)

An alternative name for `is_intersection_empty(X, Y)`.
"""
const isdisjoint = is_intersection_empty
