export is_intersection_empty,
       isdisjoint

# --- disjointness check for lazy sets ---

"""
    is_intersection_empty(X::LazySet{N},
                          Y::LazySet{N},
                          witness::Bool=false
                          ) where {N<:Real}

Check whether two sets do not intersect, and otherwise optionally compute a
witness.

### Input

- `X`       -- set
- `Y`       -- another set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ∩ Y = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``X ∩ Y = ∅``
  * `(false, v)` iff ``X ∩ Y ≠ ∅`` and ``v ∈ X ∩ Y``

### Algorithm

This is a fallback implementation that computes the concrete intersection,
`intersection`, of the given sets.

A witness is constructed using the `an_element` implementation of the result.
"""
function is_intersection_empty(X::LazySet{N},
                               Y::LazySet{N},
                               witness::Bool=false
                              ) where {N<:Real}
    cap = intersection(X, Y)
    empty_intersection = isempty(cap)
    if witness
        if empty_intersection
            return (true, N[])
        else
            return (false, an_element(cap))
        end
    end
    return empty_intersection
end

# --- alias ---

"""
    isdisjoint(X, Y)

An alternative name for `is_intersection_empty(X, Y)`.
"""
const isdisjoint = is_intersection_empty

# --- AbstractHyperrectangle ---

"""
    is_intersection_empty(H1::AbstractHyperrectangle{N},
                          H2::AbstractHyperrectangle{N},
                          witness::Bool=false
                         ) where {N<:Real}

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
                              ) where {N<:Real}
    empty_intersection = false
    center_diff = center(H2) - center(H1)
    @inbounds for i in eachindex(center_diff)
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
    @inbounds for i in eachindex(center_diff)
        if v[i] <= c2[i]
            # second center is right of first center
            v[i] += min(radius_hyperrectangle(H1, i), center_diff[i])
        else
            # second center is left of first center
            v[i] -= min(radius_hyperrectangle(H1, i), -center_diff[i])
        end
    end
    return (false, v)
end

# --- disjointness check for 1D intervals ---

"""
    is_intersection_empty(I1::Interval{N}, I2::Interval{N}, witness::Bool=false) where {N<:Real}

Check whether two intervals do not intersect, and otherwise optionally
compute a witness.

### Input

- `I1`      -- first interval
- `I2`      -- second interval
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``I1 ∩ I2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``I1 ∩ I2 = ∅``
  * `(false, v)` iff ``I1 ∩ I2 ≠ ∅`` and ``v ∈ I1 ∩ I2``

### Algorithm

``I1 ∩ I2 ≠ ∅`` iff there is a gap between the left-most point of the second
interval and the left-most point of the first interval, or vice-versa.

A witness is computed by taking the maximum over the left-most points of each
interval, which is guaranteed to belong to the intersection.
"""
function is_intersection_empty(I1::Interval{N}, I2::Interval{N}, witness::Bool=false) where {N<:Real}
    if witness
        return _is_intersection_empty(I1, I2, Val(true))
    else
        return _is_intersection_empty(I1, I2, Val(false))
    end
end

function _is_intersection_empty(I1::Interval{N}, I2::Interval{N}, witness::Val{false}) where {N<:Real}
    return !_leq(min(I2), max(I1)) || !_leq(min(I1), max(I2))
end

function _is_intersection_empty(I1::Interval{N}, I2::Interval{N}, witness::Val{true}) where {N<:Real}
    check  = _is_intersection_empty(I1, I2, Val(false))
    return (check, [max(min(I1), min(I2))])
end

# --- AbstractSingleton ---

# common code for singletons
@inline function is_intersection_empty_helper_singleton(
        S::AbstractSingleton{N}, X::LazySet{N}, witness::Bool=false
       ) where {N<:Real}
    empty_intersection = element(S) ∉ X
    if witness
        return (empty_intersection, empty_intersection ? N[] : element(S))
    else
        return empty_intersection
    end
end

"""
    is_intersection_empty(X::LazySet{N},
                          S::AbstractSingleton{N},
                          witness::Bool=false
                         ) where {N<:Real}

Check whether a convex set and a singleton do not intersect, and otherwise
optionally compute a witness.

### Input

- `X`       -- convex set
- `S`       -- singleton
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ∩ X = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``S ∩ X = ∅``
  * `(false, v)` iff ``S ∩ X ≠ ∅`` and
    `v` = `element(S)` ``∈ S ∩ X``

### Algorithm

``S ∩ X = ∅`` iff `element(S)` ``∉ X``.
"""
function is_intersection_empty(X::LazySet{N},
                               S::AbstractSingleton{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, X, witness)
end

# symmetric method
function is_intersection_empty(S::AbstractSingleton{N},
                               X::LazySet{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, X, witness)
end

"""
    is_intersection_empty(S1::AbstractSingleton{N},
                          S2::AbstractSingleton{N},
                          witness::Bool=false
                         ) where {N<:Real}

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
                              ) where {N<:Real}
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
                         ) where {N<:Real}

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
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, H, witness)
end

# symmetric method
function is_intersection_empty(S::AbstractSingleton{N},
                               H::AbstractHyperrectangle{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, H, witness)
end


# --- Ball2 ---


"""
    is_intersection_empty(B1::Ball2{N},
                          B2::Ball2{N},
                          witness::Bool=false
                         ) where {N<:AbstractFloat}

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
                              ) where {N<:AbstractFloat}
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
                         ) where {N<:Real}

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

For witness production we fall back to a less efficient implementation for
general sets as the first argument.
"""
function is_intersection_empty(Z::Zonotope{N},
                               H::Hyperplane{N},
                               witness::Bool=false
                              ) where {N<:Real}
    if witness
        # use less efficient implementation that supports witness production
        return invoke(is_intersection_empty,
                      Tuple{LazySet{N}, Hyperplane{N}, Bool},
                      Z, H, witness)
    end

    v = H.b - dot(H.a, Z.center)
    p = ngens(Z)
    if p == 0
        abs_sum = zero(N)
    else
        abs_sum = sum(abs(dot(H.a, Z.generators[:, i])) for i = 1:ngens(Z))
    end
    return v < -abs_sum || v > abs_sum
end

# symmetric method
function is_intersection_empty(H::Hyperplane{N},
                               Z::Zonotope{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty(Z, H, witness)
end

"""
    is_intersection_empty(Z1::Zonotope{N}, Z2::Zonotope{N}, witness::Bool=false
                         ) where {N<:Real}

Check whether two zonotopes do not intersect, and otherwise optionally compute a
witness.

### Input

- `Z1`      -- zonotope
- `Z2`      -- zonotope
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``Z1 ∩ Z2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``Z1 ∩ Z2 = ∅``
  * `(false, v)` iff ``Z1 ∩ Z2 ≠ ∅`` and ``v ∈ Z1 ∩ Z2``

### Algorithm

``Z1 ∩ Z2 ≠ ∅`` iff ``c_1 - c_2 ∈ Z(0, (g_1, g_2))`` where ``c_i`` and ``g_i``
are the center and generators of zonotope `Zi` and ``Z(c, g)`` represents the
zonotope with center ``c`` and generators ``g``.
"""
function is_intersection_empty(Z1::Zonotope{N}, Z2::Zonotope{N},
                               witness::Bool=false) where {N<:Real}
    n = dim(Z1)
    @assert n == dim(Z2) "zonotopes need to have the same dimensions"
    Z = Zonotope(zeros(N, n), hcat(Z1.generators, Z2.generators))
    result = (center(Z1) - center(Z2)) ∈ Z
    if result
        return witness ? (true, N[]) : true
    elseif witness
        error("witness production is not supported yet")
    else
        return false
    end
end

"""
    is_intersection_empty(ls1::LineSegment{N},
                          ls2::LineSegment{N},
                          witness::Bool=false
                         ) where {N<:Real}

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
                              ) where {N<:Real}
    function cross(ls1::Vector{N}, ls2::Vector{N}) where {N<:Real}
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


# --- Hyperplane ---


@inline function is_intersection_empty_helper_hyperplane(
        hp::Union{Hyperplane{N}, Line2D{N}},
        X::LazySet{N},
        witness::Bool=false
       ) where {N<:Real}
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


"""
    is_intersection_empty(X::LazySet{N},
                          hp::Union{Hyperplane{N}, Line2D{N}},
                          [witness]::Bool=false
                         ) where {N<:Real}

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
                               hp::Union{Hyperplane{N}, Line2D{N}},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_hyperplane(hp, X, witness)
end

# symmetric function
function is_intersection_empty(hp::Union{Hyperplane{N}, Line2D{N}},
                               X::LazySet{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_hyperplane(hp, X, witness)
end

# disambiguation
function is_intersection_empty(hp1::Union{Hyperplane{N}, Line2D{N}},
                               hp2::Union{Hyperplane{N}, Line2D{N}},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_hyperplane(hp1, hp2, witness)
end

# disambiguation
function is_intersection_empty(hp::Union{Hyperplane{N}, Line2D{N}},
                               S::AbstractSingleton{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, hp, witness)
end

# symmetric method
function is_intersection_empty(S::AbstractSingleton{N},
                               hp::Union{Hyperplane{N}, Line2D{N}},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, hp, witness)
end


# --- HalfSpace ---


@inline function is_intersection_empty_helper_halfspace(
        hs::HalfSpace{N},
        X::LazySet{N},
        witness::Bool=false
       ) where {N<:Real}
    if !witness
        return -ρ(-hs.a, X) > hs.b
    end

    # for witness production, we compute the support vector instead
    svec = σ(-hs.a, X)
    empty_intersection = svec ∉ hs
    v = empty_intersection ? N[] : svec
    return (empty_intersection, v)
end

"""
    is_intersection_empty(X::LazySet{N},
                          hs::HalfSpace{N},
                          [witness]::Bool=false
                         ) where {N<:Real}

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
`X` is a lazy intersection, options can be passed to the line search algorithm.
"""
function is_intersection_empty(X::LazySet{N},
                               hs::HalfSpace{N},
                               witness::Bool=false
                               ) where {N<:Real}
    return is_intersection_empty_helper_halfspace(hs, X, witness)
end

# symmetric function
function is_intersection_empty(hs::HalfSpace{N},
                               X::LazySet{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_halfspace(hs, X, witness)
end

"""
    is_intersection_empty(hs1::HalfSpace{N},
                          hs2::HalfSpace{N},
                          [witness]::Bool=false
                         ) where {N<:Real}

Check whether two half-spaces do not intersect, and otherwise optionally compute
a witness.

### Input

- `hs1`     -- half-space
- `hs2`     -- half-space
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``hs1 ∩ hs2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``hs1 ∩ hs2 = ∅``
  * `(false, v)` iff ``hs1 ∩ hs2 ≠ ∅`` and ``v ∈ hs1 ∩ hs2``

### Algorithm

Two half-spaces do not intersect if and only if their normal vectors point in
the opposite direction and there is a gap between the two defining hyperplanes.

The latter can be checked as follows:
Let ``hs_1 : a_1⋅x = b_1`` and ``hs2 : a_2⋅x = b_2``.
Then we already know that ``a_2 = -k⋅a_1`` for some positive scaling factor
``k``.
Let ``x_1`` be a point on the defining hyperplane of ``hs_1``.
We construct a line segment from ``x_1`` to the point ``x_2`` on the defining
hyperplane of ``hs_2`` by shooting a ray from ``x_1`` with direction ``a_1``.
Thus we look for a factor ``s`` such that ``(x_1 + s⋅a_1)⋅a_2 = b_2``.
This gives us ``s = (b_2 - x_1⋅a_2) / (-k a_1⋅a_1)``.
The gap exists if and only if ``s`` is positive.

If the normal vectors do not point in opposite directions, then the defining
hyperplanes intersect and we can produce a witness as follows.
All points ``x`` in this intersection satisfy ``a_1⋅x = b_1`` and ``a_2⋅x = b_2``.
Thus we have ``(a_1 + a_2)⋅x = b_1+b_2``.
We now find a dimension where ``a_1 + a_2`` is non-zero, say, ``i``.
Then the result is a vector with one non-zero entry in dimension ``i``, defined
as ``[0, …, 0, (b_1 + b_2)/(a_1[i] + a_2[i]), 0, …, 0]``.
Such a dimension ``i`` always exists.
"""
function is_intersection_empty(hs1::HalfSpace{N},
                               hs2::HalfSpace{N},
                               witness::Bool=false
                              ) where {N<:Real}
    a1 = hs1.a
    a2 = hs2.a
    issamedir, k = samedir(a1, -a2)
    if issamedir
        x1 = an_element(Hyperplane(a1, hs1.b))
        b2 = hs2.b
        s = (b2 - dot(x1, a2)) / (-k * dot(a1, a1))
        empty_intersection = s > 0
        if witness
            if empty_intersection
                v = N[]
            else
                # both defining hyperplanes are contained in each half-space
                v = x1
            end
        end
    else
        empty_intersection = false
        if witness
            v = zeros(N, length(a1))
            for i in 1:length(a1)
                a_sum_i = a1[i] + a2[i]
                if a_sum[i] != 0
                    v[i] = (hs1.b + hs2.b) / a_sum_i
                    break
                end
            end
        end
    end

    if !witness
        return empty_intersection
    else
        return (empty_intersection, v)
    end
end

# disambiguation
function is_intersection_empty(H::HalfSpace{N},
                               S::AbstractSingleton{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, H, witness)
end

# symmetric method
function is_intersection_empty(S::AbstractSingleton{N},
                               H::HalfSpace{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, H, witness)
end

# disambiguation
function is_intersection_empty(hp::Union{Hyperplane{N}, Line2D{N}},
                               hs::HalfSpace{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_halfspace(hs, hp, witness)
end

# symmetric method
function is_intersection_empty(hs::HalfSpace{N},
                               hp::Union{Hyperplane{N}, Line2D{N}},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_halfspace(hs, hp, witness)
end


# --- polyhedra ---


"""
    is_intersection_empty(P::AbstractPolyhedron{N},
                          X::LazySet{N},
                          witness::Bool=false;
                          solver=default_lp_solver(N)
                         ) where {N<:Real}

Check whether two polyhedra do not intersect.

### Input

- `P`         -- polyhedron
- `X`         -- another set (see the Notes section below)
- `witness`   -- (optional, default: `false`) compute a witness if activated
- `solver`    -- (optional, default: `default_lp_solver(N)`) the backend used to
                 solve the linear program
- `algorithm` -- (optional, default: `"exact"`) algorithm keyword, one of:
                 * `"exact" (exact, uses a feasibility LP)
                 * `"sufficient" (sufficient, uses half-space checks)

### Output

* If `witness` option is deactivated: `true` iff ``P ∩ X = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``P ∩ X = ∅``
  * `(false, v)` iff ``P ∩ X ≠ ∅`` and ``v ∈ P ∩ X``

### Notes

For `algorithm == "exact"`, we assume that `constraints_list(X)` is defined.
For `algorithm == "sufficient"`, witness production is not supported.

### Algorithm

For `algorithm == "exact"`, see [`isempty(P::HPoly, ::Bool)`](@ref).

For `algorithm == "sufficient"`, we rely on the intersection check between the
set `X` and each constraint in `P`.
This means one support function evaluation of `X` for each constraint of `P`.
With the sufficiency algorithm, this function may return `false` even in the case
where the intersection is empty. On the other hand, if the algorithm returns
`true`, then it is guaranteed that the intersection is empty.
"""
function is_intersection_empty(P::AbstractPolyhedron{N},
                               X::LazySet{N},
                               witness::Bool=false;
                               solver=default_lp_solver(N),
                               algorithm="exact"
                              ) where {N<:Real}
    if algorithm == "sufficient"
        # sufficient check for empty intersection using half-space checks
        for Hi in constraints_list(P)
            if is_intersection_empty(X, Hi)
                if witness
                    return (true, N[])
                end
                return true
            end
        end
        if witness
            error("witness production is not supported yet")
        end
        return false
    elseif algorithm == "exact"
        # exact check for empty intersection using a feasibility LP
        clist_P = _normal_Vector(P) # TODO
        clist_X =_normal_Vector(X) # TODO
        return isempty(HPolyhedron([clist_P; clist_X]), witness; solver=solver)
    else
        error("algorithm $algorithm unknown")
    end
end

# symmetric method
function is_intersection_empty(X::LazySet{N},
                               P::AbstractPolyhedron{N},
                               witness::Bool=false;
                               solver=default_lp_solver(N),
                               algorithm="exact"
                              ) where {N<:Real}
    return is_intersection_empty(P, X, witness;
                                 solver=solver, algorithm=algorithm)
end

# disambiguation
function is_intersection_empty(P::AbstractPolyhedron{N},
                               Q::AbstractPolyhedron{N},
                               witness::Bool=false;
                               solver=default_lp_solver(N),
                               algorithm="exact"
                              ) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{typeof(P), LazySet{N}, Bool},
                  P, Q, witness; solver=solver, algorithm=algorithm)
end

# disambiguation
function is_intersection_empty(P::AbstractPolyhedron{N},
                               hs::HalfSpace{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_halfspace(hs, P, witness)
end

# symmetric method
function is_intersection_empty(hs::HalfSpace{N},
                               P::AbstractPolyhedron{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_halfspace(hs, P, witness)
end

# disambiguation
function is_intersection_empty(P::AbstractPolyhedron{N},
                               S::AbstractSingleton{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, P, witness)
end

# symmetric method
function is_intersection_empty(S::AbstractSingleton{N},
                               P::AbstractPolyhedron{N},
                               witness::Bool=false
                              ) where {N<:Real}
    return is_intersection_empty_helper_singleton(S, P, witness)
end

# disambiguation
function is_intersection_empty(P::AbstractPolyhedron{N},
                               hp::Union{Hyperplane{N}, Line2D{N}},
                               witness::Bool=false;
                               solver=default_lp_solver(N),
                               algorithm="exact"
                              ) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{typeof(P), LazySet{N}, Bool},
                  P, hp, witness, solver=solver, algorithm=algorithm)
end

# symmetric method
function is_intersection_empty(hp::Union{Hyperplane{N}, Line2D{N}},
                               P::AbstractPolyhedron{N},
                               witness::Bool=false;
                               solver=default_lp_solver(N),
                               algorithm="exact"
                              ) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{typeof(P), LazySet{N}, Bool},
                  P, hp, witness, solver=solver, algorithm=algorithm)
end


# --- union ---


"""
    is_intersection_empty(cup::UnionSet{N}, X::LazySet{N},
                          [witness]::Bool=false) where {N<:Real}

Check whether a union of two convex sets and another set do not intersect.

### Input

- `cup` -- union of two convex sets
- `X`   -- another set

### Output

`true` iff ``\\text{cup} ∩ X = ∅``.
"""
function is_intersection_empty(cup::UnionSet{N}, X::LazySet{N},
                               witness::Bool=false) where {N<:Real}
    return is_intersection_empty(UnionSetArray([cup.X, cup.Y]), X, witness)
end

# symmetric method
function is_intersection_empty(X::LazySet{N}, cup::UnionSet{N},
                               witness::Bool=false) where {N<:Real}
    return is_intersection_empty(cup, X, witness)
end

# disambiguation
function is_intersection_empty(cup1::UnionSet{N}, cup2::UnionSet{N},
                               witness::Bool=false) where {N<:Real}
    return is_intersection_empty(UnionSetArray([cup1.X, cup1.Y]),
                                 UnionSetArray([cup2.X, cup2.Y]), witness)
end

"""
    is_intersection_empty(cup::UnionSetArray{N}, X::LazySet{N},
                          [witness]::Bool=false) where {N<:Real}

Check whether a union of a finite number of convex sets and another set do not
intersect.

### Input

- `cup` -- union of a finite number of convex sets
- `X`   -- another set

### Output

`true` iff ``\\text{cup} ∩ X = ∅``.
"""
function is_intersection_empty(cup::UnionSetArray{N}, X::LazySet{N},
                               witness::Bool=false) where {N<:Real}
    result = true
    w = N[]
    for Y in array(cup)
        if witness
            result, w = is_intersection_empty(Y, X, witness)
        else
            result = is_intersection_empty(Y, X, witness)
        end
        if !result
            break
        end
    end
    return witness ? (result, w) : result
end

# symmetric method
function is_intersection_empty(X::LazySet{N}, cup::UnionSetArray{N},
                               witness::Bool=false) where {N<:Real}
    return is_intersection_empty(cup, X, witness)
end

# disambiguation
function is_intersection_empty(cup1::UnionSet{N}, cup2::UnionSetArray{N},
                               witness::Bool=false) where {N<:Real}
    return is_intersection_empty(UnionSetArray([cup1.X, cup1.Y]), cup2, witness)
end

# disambiguation
function is_intersection_empty(cup1::UnionSetArray{N}, cup2::UnionSet{N},
                               witness::Bool=false) where {N<:Real}
    return is_intersection_empty(cup1, UnionSetArray([cup2.X, cup2.Y]), witness)
end

# disambiguation
function is_intersection_empty(cup1::UnionSetArray{N},
                               cup2::UnionSetArray{N},
                               witness::Bool=false) where {N<:Real}
    result = true
    w = N[]
    for X in array(cup1)
        for Y in array(cup2)
            if witness
                result, w = is_intersection_empty(Y, X, witness)
            else
                result = is_intersection_empty(Y, X, witness)
            end
            if !result
                break
            end
        end
    end
    return witness ? (result, w) : result
end


# --- Universe ---


"""
    is_intersection_empty(U::Universe{N}, X::LazySet{N},
                          [witness]::Bool=false) where {N<:Real}

Check whether a universe and another set do not intersect.

### Input

- `U` -- universe
- `X` -- another set

### Output

`true` iff ``X ≠ ∅``.
"""
function is_intersection_empty(U::Universe{N}, X::LazySet{N}, witness::Bool=false) where {N<:Real}
    return _is_intersection_empty_universe(X, U, witness)
end

# symmetric method
function is_intersection_empty(X::LazySet{N}, U::Universe{N}, witness::Bool=false) where {N<:Real}
    return _is_intersection_empty_universe(X, U, witness)
end

function _is_intersection_empty_universe(X, U::Universe{N}, witness) where {N}
    @assert dim(X) == dim(U) "the dimensions of the given sets should match but they are $(dim(X)) and $(dim(U)) respectively"
    result = isempty(X)
    if result
        return witness ? (result, N[]) : result
    else
        return witness ? (result, an_element(X)) : result
    end
end

# disambiguation
function is_intersection_empty(U::Universe{N}, ::Universe{N},
                               witness::Bool=false) where {N<:Real}
    return witness ? (false, an_element(U)) : false
end
function is_intersection_empty(P::AbstractPolyhedron{N}, U::Universe{N},
                               witness::Bool=false) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{Universe{N}, LazySet{N}, Bool},
                  U, P, witness)
end
function is_intersection_empty(U::Universe{N}, P::AbstractPolyhedron{N},
                               witness::Bool=false) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{Universe{N}, LazySet{N}, Bool},
                  U, P, witness)
end
function is_intersection_empty(S::AbstractSingleton{N}, U::Universe{N},
                               witness::Bool=false) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{Universe{N}, LazySet{N}, Bool},
                  U, S, witness)
end
function is_intersection_empty(U::Universe{N}, S::AbstractSingleton{N},
                               witness::Bool=false) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{Universe{N}, LazySet{N}, Bool},
                  U, S, witness)
end
function is_intersection_empty(hs::HalfSpace{N}, U::Universe{N},
                               witness::Bool=false) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{Universe{N}, LazySet{N}, Bool},
                  U, hs, witness)
end
function is_intersection_empty(U::Universe{N},
                               hp::Union{Hyperplane{N}, Line2D{N}},
                               witness::Bool=false) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{Universe{N}, LazySet{N}, Bool},
                  U, hp, witness)
end
function is_intersection_empty(hp::Union{Hyperplane{N}, Line2D{N}},
                               U::Universe{N},
                               witness::Bool=false) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{Universe{N}, LazySet{N}, Bool},
                  U, hp, witness)
end
function is_intersection_empty(U::Universe{N}, hs::HalfSpace{N},
                               witness::Bool=false) where {N<:Real}
    return invoke(is_intersection_empty,
                  Tuple{Universe{N}, LazySet{N}, Bool},
                  U, hs, witness)
end


# --- Complement ---


"""
    is_intersection_empty(C::Complement{N},
                          X::LazySet{N},
                          [witness]::Bool=false
                         ) where {N<:Real}

Check whether the complement of a convex set and another set do not intersect.

### Input

- `C` -- complement of a convex set
- `X` -- convex set

### Output

* If `witness` option is deactivated: `true` iff ``X ∩ C = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``X ∩ C = ∅``
  * `(false, v)` iff ``X ∩ C ≠ ∅`` and ``v ∈ X ∩ C``

### Algorithm

We fall back to `X ⊆ C.X`, which can be justified as follows:

```math
    X ∩ Y^C = ∅ ⟺ X ⊆ Y
```
"""
function is_intersection_empty(C::Complement{N}, X::LazySet{N},
                               witness::Bool=false) where {N<:Real}
    return ⊆(X, C.X, witness)
end

# symmetric method
function is_intersection_empty(X::LazySet{N}, C::Complement{N},
                               witness::Bool=false) where {N<:Real}
    return is_intersection_empty(C, X, witness)
end

"""
    is_intersection_empty(cpa::CartesianProductArray{N},
                          P::AbstractPolyhedron{N}) where {N<:Real}

Check whether a polytopic Cartesian product array intersects with a polyhedron.

### Input

- `cpa` -- Cartesian product array of polytopes
- `P`   -- polyhedron

### Output

`true` iff ``\\text{cpa} ∩ Y = ∅ ``.

### Algorithm

We first identify the blocks of `cpa` in which `P` is constrained.
Then we project `cpa` to those blocks and convert the result to an `HPolytope`
`Q`.
Finally we determine whether `Q` and the projected `P` intersect.
"""
function is_intersection_empty(cpa::CartesianProductArray{N},
                               P::AbstractPolyhedron{N}) where {N<:Real}
    cpa_low_dim, vars, _block_structure = get_constrained_lowdimset(cpa, P)
    hpoly_low_dim = HPolytope(constraints_list(cpa_low_dim))
    return isdisjoint(hpoly_low_dim, project(P, vars))
end

# symmetric method
function is_intersection_empty(P::AbstractPolyhedron{N},
                               cpa::CartesianProductArray{N}) where {N<:Real}
    return is_intersection_empty(cpa, P)
end

# disambiguation
function is_intersection_empty(cpa::CartesianProductArray{N},
                               hs::HalfSpace{N}) where {N<:Real}
    return is_intersection_empty_helper_halfspace(hs, cpa)
end
function is_intersection_empty(hs::HalfSpace{N},
                               cpa::CartesianProductArray{N}) where {N<:Real}
    return is_intersection_empty(cpa, hs)
end

function is_intersection_empty(cpa::CartesianProductArray{N,S}, U::Universe{N},
                               witness::Bool=false) where {N<:Real, S<:LazySet{N}}
    return _is_intersection_empty_universe(cpa, U, witness)
end
function is_intersection_empty(U::Universe{N}, cpa::CartesianProductArray{N, S},
                               witness::Bool=false) where {N<:Real, S<:LazySet{N}}
    return _is_intersection_empty_universe(cpa, U, witness)
end

"""
    is_intersection_empty(X::CartesianProductArray{N},
                          Y::CartesianProductArray{N}) where {N<:Real}

Check whether two Cartesian products of a finite number of convex sets do not
intersect.

### Input

- `X` -- Cartesian product array of convex sets
- `Y` -- Cartesian product array of convex sets

### Output

`true` iff ``X ∩ Y = ∅``.
"""
function is_intersection_empty(X::CartesianProductArray{N},
                               Y::CartesianProductArray{N}) where {N<:Real}
    @assert same_block_structure(array(X), array(Y)) "block structure has to " *
        "be the same"

    for i in 1:length(X.array)
        if isdisjoint(X.array[i], Y.array[i])
            return true
        end
    end
    return false
end

"""
    is_intersection_empty(cpa::CartesianProductArray{N},
                          H::AbstractHyperrectangle{N}) where {N<:Real}

Check whether a Cartesian product of a finite number of convex sets and a
hyperrectangular set do not intersect, and otherwise optionally compute a
witness.

### Input

- `cpa`     -- Cartesian product of a finite number of convex sets
- `H`       -- hyperrectangular set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``cpa ∩ H = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``cpa ∩ H = ∅``
  * `(false, v)` iff ``cpa ∩ H ≠ ∅`` and ``v ∈ cpa ∩ H``

### Algorithm

The sets `cpa` and `H` are disjoint if and only if at least one block of `cpa`
and the corresponding projection of `H` are disjoint.
We perform these checks sequentially.
"""
function is_intersection_empty(cpa::CartesianProductArray{N},
                               H::AbstractHyperrectangle{N},
                               witness::Bool=false) where {N<:Real}
    if witness
        w = zeros(N, dim(H))
    end
    n = dim(H)
    block_start = 1
    for bi in array(cpa)
        ni = dim(bi)
        block = block_start:(block_start + ni - 1)
        Hi = project(H, block, Hyperrectangle, n)
        res = isdisjoint(bi, Hi, witness)
        if witness
            if res[1]
                return (true, N[])
            else
                w[block] = res[2]
            end
        elseif res
            return true
        end
        block_start += ni
    end
    return witness ? (false, w) : false
end

# symmetric method
function is_intersection_empty(H::AbstractHyperrectangle{N},
                               cpa::CartesianProductArray{N},
                               witness::Bool=false) where {N<:Real}
    return is_intersection_empty(cpa, H, witness)
end
