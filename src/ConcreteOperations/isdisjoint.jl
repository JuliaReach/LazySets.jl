"""
    isdisjoint(X::LazySet, Y::LazySet, [witness]::Bool=false)

Check whether two sets do not intersect, and otherwise optionally compute a
witness.

### Input

- `X`       -- set
- `Y`       -- set
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
function isdisjoint(X::LazySet, Y::LazySet, witness::Bool=false)
    if _isdisjoint_convex_sufficient(X, Y)
        return _witness_result_empty(witness, true, X, Y)
    end
    return _isdisjoint_general(X, Y, witness)
end

function _isdisjoint_general(X::LazySet, Y::LazySet, witness::Bool=false)
    cap = intersection(X, Y)
    empty_intersection = isempty(cap)
    if empty_intersection
        return _witness_result_empty(witness, true, X, Y)
    end
    return witness ? (false, an_element(cap)) : false
end

# quick sufficient check that tries to find a separating hyperplane
# the result `true` is also sufficient for non-convex sets
function _isdisjoint_convex_sufficient(X::LazySet, Y::LazySet)
    x = an_element(X)
    y = an_element(Y)
    d = x - y
    return ρ(d, Y) < -ρ(-d, X)
end

# conversion for IA types
isdisjoint(X::LazySet, Y::IA.Interval, witness::Bool=false) = isdisjoint(X, Interval(Y), witness)
isdisjoint(X::IA.Interval, Y::LazySet, witness::Bool=false) = isdisjoint(Interval(X), Y, witness)

function isdisjoint(X::LazySet, Y::IA.IntervalBox, witness::Bool=false)
    return isdisjoint(X, convert(Hyperrectangle, Y), witness)
end
function isdisjoint(X::IA.IntervalBox, Y::LazySet, witness::Bool=false)
    return isdisjoint(convert(Hyperrectangle, X), Y, witness)
end

"""
    isdisjoint(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle,
               [witness]::Bool=false)

Check whether two hyperrectangular sets do not intersect, and otherwise
optionally compute a witness.

### Input

- `H1`      -- hyperrectangular set
- `H2`      -- hyperrectangular set
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
function isdisjoint(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle,
                    witness::Bool=false)
    empty_intersection = false
    center_diff = center(H2) - center(H1)
    @inbounds for i in eachindex(center_diff)
        if abs(center_diff[i]) >
           radius_hyperrectangle(H1, i) + radius_hyperrectangle(H2, i)
            empty_intersection = true
            break
        end
    end

    if empty_intersection
        return _witness_result_empty(witness, true, H1, H2)
    elseif !witness
        return false
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

# common code for singletons
function _isdisjoint_singleton(S::AbstractSingleton, X::LazySet,
                               witness::Bool=false)
    s = element(S)
    empty_intersection = s ∉ X
    return _witness_result_empty(witness, empty_intersection, S, X, s)
end

"""
    isdisjoint(X::LazySet, S::AbstractSingleton, [witness]::Bool=false)

Check whether a set and a set with a single value do not intersect, and
otherwise optionally compute a witness.

### Input

- `X`       -- set
- `S`       -- set with a single value
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S ∩ X = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``S ∩ X = ∅``
  * `(false, v)` iff ``S ∩ X ≠ ∅`` and `v` = `element(S)` ``∈ S ∩ X``

### Algorithm

``S ∩ X = ∅`` iff `element(S)` ``∉ X``.
"""
@commutative function isdisjoint(X::LazySet, S::AbstractSingleton,
                                 witness::Bool=false)
    return _isdisjoint_singleton(S, X, witness)
end

# disambiguations
for ST in [:AbstractPolyhedron, :AbstractZonotope, :AbstractHyperrectangle,
           :Hyperplane, :Line2D, :HalfSpace, :CartesianProductArray, :UnionSet,
           :UnionSetArray]
    @eval @commutative function isdisjoint(X::($ST), S::AbstractSingleton,
                                           witness::Bool=false)
        return _isdisjoint_singleton(S, X, witness)
    end
end

"""
    isdisjoint(S1::AbstractSingleton, S2::AbstractSingleton,
               [witness]::Bool=false)

Check whether two sets with a single value do not intersect, and otherwise
optionally compute a witness.

### Input

- `S1` -- set with a single value
- `S2` -- set with a single value
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``S1 ∩ S2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``S1 ∩ S2 = ∅``
  * `(false, v)` iff ``S1 ∩ S2 ≠ ∅`` and `v` = `element(S1)` ``∈ S1 ∩ S2``

### Algorithm

``S1 ∩ S2 = ∅`` iff ``S1 ≠ S2``.
"""
function isdisjoint(S1::AbstractSingleton, S2::AbstractSingleton,
                    witness::Bool=false)
    s1 = element(S1)
    empty_intersection = !isapprox(s1, element(S2))
    return _witness_result_empty(witness, empty_intersection, S1, S2, s1)
end

"""
    isdisjoint(B1::Ball2, B2::Ball2, [witness]::Bool=false)

Check whether two balls in the 2-norm do not intersect, and otherwise optionally
compute a witness.

### Input

- `B1`      -- ball in the 2-norm
- `B2`      -- ball in the 2-norm
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
function isdisjoint(B1::Ball2, B2::Ball2, witness::Bool=false)
    center_diff_normed = norm(center(B2) - center(B1), 2)
    empty_intersection = center_diff_normed > B1.radius + B2.radius

    if empty_intersection
        return _witness_result_empty(witness, true, B1, B2)
    elseif !witness
        return false
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

"""
    isdisjoint(Z::AbstractZonotope, H::Hyperplane, [witness]::Bool=false)

Check whether a zonotopic set and a hyperplane do not intersect, and otherwise
optionally compute a witness.

### Input

- `Z`       -- zonotopic set
- `H`       -- hyperplane
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
@commutative function isdisjoint(Z::AbstractZonotope, H::Hyperplane,
                                 witness::Bool=false)
    return _isdisjoint_zonotope_hyperplane(Z, H, witness)
end

@commutative function isdisjoint(Z::AbstractZonotope, H::Line2D,
                                 witness::Bool=false)
    return _isdisjoint_zonotope_hyperplane(Z, H, witness)
end

function _isdisjoint_zonotope_hyperplane(Z::AbstractZonotope,
                                         H::Union{Hyperplane,Line2D},
                                         witness::Bool=false)
    if witness
        return _isdisjoint(Z, H, Val(true))
    else
        return _isdisjoint(Z, H, Val(false))
    end
end

function _isdisjoint(Z::AbstractZonotope, H::Union{Hyperplane,Line2D},
                     ::Val{false})
    c, G = center(Z), genmat(Z)
    v = H.b - dot(H.a, c)

    p = size(G, 2)
    p == 0 && return !isapproxzero(v)
    asum = abs_sum(H.a, G)
    return !_geq(v, -asum) || !_leq(v, asum)
end

function _isdisjoint(Z::AbstractZonotope, H::Union{Hyperplane,Line2D},
                     ::Val{true})
    return _isdisjoint_hyperplane(H, Z, true)
end

"""
    isdisjoint(Z1::AbstractZonotope, Z2::AbstractZonotope,
               [witness]::Bool=false; [solver]=nothing)

Check whether two zonotopic sets do not intersect, and otherwise optionally
compute a witness.

### Input

- `Z1`      -- zonotopic set
- `Z2`      -- zonotopic set
- `witness` -- (optional, default: `false`) compute a witness if activated
- `solver`  -- (optional, default: `nothing`) the backend used to solve the
               linear program

### Output

* If `witness` option is deactivated: `true` iff ``Z1 ∩ Z2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``Z1 ∩ Z2 = ∅``
  * `(false, v)` iff ``Z1 ∩ Z2 ≠ ∅`` and ``v ∈ Z1 ∩ Z2``

### Algorithm

The algorithm is taken from [1].

``Z1 ∩ Z2 = ∅`` iff ``c_1 - c_2 ∉ Z(0, (g_1, g_2))`` where ``c_i`` and ``g_i``
are the center and generators of zonotope `Zi` and ``Z(c, g)`` represents the
zonotope with center ``c`` and generators ``g``.

[1] L. J. Guibas, A. T. Nguyen, L. Zhang: *Zonotopes as bounding volumes*. SODA
2003.
"""
function isdisjoint(Z1::AbstractZonotope, Z2::AbstractZonotope,
                    witness::Bool=false; solver=nothing)
    if _isdisjoint_convex_sufficient(Z1, Z2)
        return _witness_result_empty(witness, true, Z1, Z2)
    end

    n = dim(Z1)
    @assert n == dim(Z2) "the zonotopes need to have the same dimensions"
    N = promote_type(eltype(Z1), eltype(Z2))
    Z = Zonotope(zeros(N, n), hcat(genmat(Z1), genmat(Z2)))
    result = !∈(center(Z1) - center(Z2), Z; solver=solver)
    if result
        return _witness_result_empty(witness, true, N)
    elseif witness
        error("witness production is not supported yet")
    else
        return false
    end
end

"""
    isdisjoint(L1::LineSegment, L2::LineSegment, [witness]::Bool=false)

Check whether two line segments do not intersect, and otherwise optionally
compute a witness.

### Input

- `L1`      -- line segment
- `L2`      -- line segment
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``L1 ∩ L2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``L1 ∩ L2 = ∅``
  * `(false, v)` iff ``L1 ∩ L2 ≠ ∅`` and ``v ∈ L1 ∩ L2``

### Algorithm

The algorithm is inspired from [here](https://stackoverflow.com/a/565282), which
again is the special 2D case of a 3D algorithm from [1].

We first check if the two line segments are parallel, and if so, if they are
collinear. In the latter case, we check membership of any of the end points in
the other line segment. Otherwise the lines are not parallel, so we can solve an
equation of the intersection point, if it exists.

[1] Ronald Goldman. *Intersection of two lines in three-space*. Graphics Gems
1990.
"""
function isdisjoint(L1::LineSegment, L2::LineSegment, witness::Bool=false)
    r = L1.q - L1.p
    if all(isapproxzero, r)
        # first line segment is a point
        empty_intersection = L1.q ∉ L2
        return _witness_result_empty(witness, empty_intersection, L1, L2, L1.q)
    end

    s = L2.q - L2.p
    if all(isapproxzero, s)
        # second line segment is a point
        empty_intersection = L2.q ∉ L1
        return _witness_result_empty(witness, empty_intersection, L1, L2, L2.q)
    end

    p1p2 = L2.p - L1.p
    u_numerator = right_turn(p1p2, r)
    u_denominator = right_turn(r, s)

    if u_denominator == 0
        # line segments are parallel
        if u_numerator == 0
            # line segments are collinear
            if L1.p ∈ L2
                empty_intersection = false
                if witness
                    v = L1.p
                end
            elseif L1.q ∈ L2
                empty_intersection = false
                if witness
                    v = L1.q
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
            t = right_turn(p1p2, s) / u_denominator
            empty_intersection = t < 0 || t > 1
            if witness
                v = L1.p + t * r
            end
        end
    end
    if witness && !empty_intersection
        return (false, v)
    end
    return _witness_result_empty(witness, empty_intersection, L1, L2)
end

function _isdisjoint_hyperplane(hp::Union{Hyperplane,Line2D}, X::LazySet,
                                witness::Bool=false)
    if !isconvextype(typeof(X))
        error("this implementation requires a convex set")
    end

    normal_hp = hp.a
    sv_left = σ(-normal_hp, X)
    if -dot(sv_left, -normal_hp) <= hp.b
        sv_right = σ(normal_hp, X)
        empty_intersection = (hp.b > dot(sv_right, normal_hp))
    else
        empty_intersection = true
    end
    if !witness || empty_intersection
        return _witness_result_empty(witness, empty_intersection, hp, X)
    end
    # compute witness
    point_hp = an_element(hp)
    point_line = sv_left
    dir_line = sv_right - sv_left
    d = dot((point_hp - point_line), normal_hp) /
        dot(dir_line, normal_hp)
    v = d * dir_line + point_line
    return (false, v)
end

"""
    isdisjoint(X::LazySet, hp::Hyperplane, [witness]::Bool=false)

Check whether a convex set an a hyperplane do not intersect, and otherwise
optionally compute a witness.

### Input

- `X`       -- convex set
- `hp`      -- hyperplane
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ∩ hp = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``X ∩ hp = ∅``
  * `(false, v)` iff ``X ∩ hp ≠ ∅`` and ``v ∈ X ∩ hp``

### Algorithm

A convex set intersects with a hyperplane iff the support function in the
negative resp. positive direction of the hyperplane's normal vector ``a`` is to
the left resp. right of the hyperplane's constraint ``b``:

```math
-ρ(-a, X) ≤ b ≤ ρ(a, X)
```

For witness generation, we compute a line connecting the support vectors to the
left and right, and then take the intersection of the line with the hyperplane.
We follow
[this algorithm](https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection#Algebraic_form)
for the line-hyperplane intersection.
"""
@commutative function isdisjoint(X::LazySet, hp::Hyperplane,
                                 witness::Bool=false)
    return _isdisjoint_hyperplane(hp, X, witness)
end

@commutative function isdisjoint(X::LazySet, L::Line2D, witness::Bool=false)
    return _isdisjoint_hyperplane(L, X, witness)
end

#disambiguations
for ST in [:AbstractPolyhedron]
    @eval @commutative function isdisjoint(X::($ST), H::Hyperplane,
                                           witness::Bool=false)
        return _isdisjoint_hyperplane(H, X, witness)
    end

    @eval @commutative function isdisjoint(X::($ST), L::Line2D, witness::Bool=false)
        return _isdisjoint_hyperplane(L, X, witness)
    end
end

function isdisjoint(hp1::Hyperplane, hp2::Hyperplane, witness::Bool=false)
    return _isdisjoint_hyperplane_hyperplane(hp1, hp2, witness)
end

@commutative function isdisjoint(hp::Hyperplane, L::Line2D, witness::Bool=false)
    return _isdisjoint_hyperplane_hyperplane(hp, L, witness)
end

function _isdisjoint_hyperplane_hyperplane(hp1::Union{Hyperplane,Line2D},
                                           hp2::Union{Hyperplane,Line2D},
                                           witness::Bool=false)
    if isequivalent(hp1, hp2)
        res = false
        if witness
            w = an_element(hp1)
        end
    else
        cap = intersection(hp1, hp2)
        res = cap isa EmptySet
        if !res && witness
            w = an_element(cap)
        end
    end
    if res
        return _witness_result_empty(witness, true, hp1, hp2)
    end
    return witness ? (false, w) : false
end

function _isdisjoint_halfspace(hs::HalfSpace, X::LazySet, witness::Bool=false)
    if !witness
        return !_leq(-ρ(-hs.a, X), hs.b)
    end

    # for witness production, we compute the support vector instead
    svec = σ(-hs.a, X)
    empty_intersection = svec ∉ hs
    return _witness_result_empty(witness, empty_intersection, hs, X, svec)
end

"""
    isdisjoint(X::LazySet, hs::HalfSpace, [witness]::Bool=false)

Check whether a set an a half-space do not intersect, and otherwise optionally
compute a witness.

### Input

- `X`       -- set
- `hs`      -- half-space
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ∩ hs = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``X ∩ hs = ∅``
  * `(false, v)` iff ``X ∩ hs ≠ ∅`` and ``v ∈ X ∩ hs``

### Algorithm

A set intersects with a half-space iff the support function in the negative
direction of the half-space's normal vector ``a`` is less than the constraint
``b`` of the half-space: ``-ρ(-a, X) ≤ b``.

For compact set `X`, we equivalently have that the support vector in the
negative direction ``-a`` is contained in the half-space: ``σ(-a) ∈ hs``.
The support vector is thus also a witness if the sets are not disjoint.
"""
@commutative function isdisjoint(X::LazySet, hs::HalfSpace, witness::Bool=false)
    return _isdisjoint_halfspace(hs, X, witness)
end

"""
    isdisjoint(H1::HalfSpace, H2::HalfSpace, [witness]::Bool=false)

Check whether two half-spaces do not intersect, and otherwise optionally compute
a witness.

### Input

- `H1`     -- half-space
- `H2`     -- half-space
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``H1 ∩ H2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``H1 ∩ H2 = ∅``
  * `(false, v)` iff ``H1 ∩ H2 ≠ ∅`` and ``v ∈ H1 ∩ H2``

### Algorithm

Two half-spaces do not intersect if and only if their normal vectors point in
the opposite direction and there is a gap between the two defining hyperplanes.

The latter can be checked as follows:
Let ``H1 : a_1⋅x = b_1`` and ``H2 : a_2⋅x = b_2``.
Then we already know that ``a_2 = -k⋅a_1`` for some positive scaling factor
``k``.
Let ``x_1`` be a point on the defining hyperplane of ``H1``.
We construct a line segment from ``x_1`` to the point ``x_2`` on the defining
hyperplane of ``hs_2`` by shooting a ray from ``x_1`` with direction ``a_1``.
Thus we look for a factor ``s`` such that ``(x_1 + s⋅a_1)⋅a_2 = b_2``.
This gives us ``s = (b_2 - x_1⋅a_2) / (-k a_1⋅a_1)``.
The gap exists if and only if ``s`` is positive.

If the normal vectors do not point in opposite directions, then the defining
hyperplanes intersect and we can produce a witness as follows.
All points ``x`` in this intersection satisfy ``a_1⋅x = b_1`` and
``a_2⋅x = b_2``. Thus we have ``(a_1 + a_2)⋅x = b_1+b_2``.
We now find a dimension where ``a_1 + a_2`` is non-zero, say, ``i``.
Then the result is a vector with one non-zero entry in dimension ``i``, defined
as ``[0, …, 0, (b_1 + b_2)/(a_1[i] + a_2[i]), 0, …, 0]``.
Such a dimension ``i`` always exists.
"""
function isdisjoint(H1::HalfSpace, H2::HalfSpace, witness::Bool=false)
    a1 = H1.a
    a2 = H2.a
    N = promote_type(eltype(H1), eltype(H2))
    issamedir, k = samedir(a1, -a2)
    if issamedir
        x1 = an_element(Hyperplane(a1, H1.b))
        b2 = H2.b
        s = (b2 - dot(x1, a2)) / (-k * dot(a1, a1))
        empty_intersection = s > 0
        # if `!empty_intersection`, x1 is a witness because both defining
        # hyperplanes are contained in each half-space
        return _witness_result_empty(witness, empty_intersection, H1, H2, x1)
    elseif !witness
        return false
    end
    # compute witness
    v = zeros(N, length(a1))
    for i in eachindex(a1)
        a_sum_i = a1[i] + a2[i]
        if a_sum[i] != 0
            v[i] = (H1.b + H2.b) / a_sum_i
            break
        end
    end
    return (false, v)
end

# disambiguations
for ST in [:AbstractPolyhedron, :AbstractZonotope, :Hyperplane, :Line2D,
           :CartesianProductArray]
    @eval @commutative function isdisjoint(X::($ST), H::HalfSpace, witness::Bool=false)
        return _isdisjoint_halfspace(H, X, witness)
    end
end

"""
    isdisjoint(P::AbstractPolyhedron, X::LazySet, [witness]::Bool=false;
               [solver]=nothing, [algorithm]="exact")

Check whether a polyhedral set and another set do not intersect, and otherwise
optionally compute a witness.

### Input

- `P`         -- polyhedral set
- `X`         -- set (see the Notes section below)
- `witness`   -- (optional, default: `false`) compute a witness if activated
- `solver`    -- (optional, default: `nothing`) the backend used to solve the
                 linear program
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

For `solver == nothing`, we fall back to `default_lp_solver(N)`.

### Algorithm

For `algorithm == "exact"`, see [`isempty(P::HPoly, ::Bool)`](@ref).

For `algorithm == "sufficient"`, we rely on the intersection check between the
set `X` and each constraint in `P`.
This requires one support-function evaluation of `X` for each constraint of `P`.
With this algorithm, the method may return `false` even in the case where the
intersection is empty. On the other hand, if the algorithm returns `true`, then
it is guaranteed that the intersection is empty.
"""
@commutative function isdisjoint(P::AbstractPolyhedron, X::LazySet,
                                 witness::Bool=false; solver=nothing,
                                 algorithm="exact")
    return _isdisjoint_polyhedron(P, X, witness; solver=solver,
                                  algorithm=algorithm)
end

function _isdisjoint_polyhedron(P::AbstractPolyhedron, X::LazySet,
                                witness::Bool=false; solver=nothing,
                                algorithm="exact")
    N = promote_type(eltype(P), eltype(X))
    if algorithm == "sufficient"
        # sufficient check for empty intersection using half-space checks
        for Hi in constraints_list(P)
            if isdisjoint(X, Hi)
                return _witness_result_empty(witness, true, N)
            end
        end
        if witness
            error("witness production is not supported yet")
        end
        return false
    elseif algorithm == "exact"
        # exact check for empty intersection using a feasibility LP
        if !is_polyhedral(X)
            error("this algorithm requires a polyhedral input")
        end
        clist_P = _normal_Vector(P) # TODO
        clist_X = _normal_Vector(X) # TODO
        if isnothing(solver)
            solver = default_lp_solver(N)
        end
        return _isempty_polyhedron_lp([clist_P; clist_X], witness; solver=solver)
    else
        error("algorithm $algorithm unknown")
    end
end

# disambiguation
function isdisjoint(X::AbstractPolyhedron, P::AbstractPolyhedron, witness::Bool=false;
                    solver=nothing, algorithm="exact")
    return _isdisjoint_polyhedron(P, X, witness)
end

"""
    isdisjoint(U::UnionSet, X::LazySet, [witness]::Bool=false)

Check whether a union of two sets and another set do not intersect, and
otherwise optionally compute a witness.

### Input

- `U`       -- union of two sets
- `X`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

`true` iff ``\\text{U} ∩ X = ∅``.
"""
@commutative function isdisjoint(U::UnionSet, X::LazySet, witness::Bool=false)
    return _isdisjoint_union(U, X, witness)
end

"""
    isdisjoint(U::UnionSetArray, X::LazySet, [witness]::Bool=false)

Check whether a union of a finite number of sets and another set do not
intersect, and otherwise optionally compute a witness.

### Input

- `U`       -- union of a finite number of sets
- `X`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

`true` iff ``\\text{U} ∩ X = ∅``.
"""
@commutative function isdisjoint(U::UnionSetArray, X::LazySet,
                                 witness::Bool=false)
    return _isdisjoint_union(U, X, witness)
end

function _isdisjoint_union(cup::Union{UnionSet,UnionSetArray}, X::LazySet{N},
                           witness::Bool=false) where {N}
    for Y in cup
        if witness
            result, w = isdisjoint(Y, X, witness)
        else
            result = isdisjoint(Y, X, witness)
        end
        if !result
            return witness ? (false, w) : false
        end
    end
    return _witness_result_empty(witness, true, N)
end

# disambiguations
for ST in [:AbstractPolyhedron, :Hyperplane, :Line2D, :HalfSpace]
    @eval @commutative function isdisjoint(U::UnionSet, X::($ST), witness::Bool=false)
        return _isdisjoint_union(U, X, witness)
    end
    @eval @commutative function isdisjoint(U::UnionSetArray, X::($ST),
                                           witness::Bool=false)
        return _isdisjoint_union(U, X, witness)
    end
end

@commutative function isdisjoint(U1::UnionSet, U2::UnionSetArray,
                                 witness::Bool=false)
    return _isdisjoint_union(U1, U2, witness)
end

function isdisjoint(U1::UnionSet, U2::UnionSet, witness::Bool=false)
    return _isdisjoint_union(U1, U2, witness)
end

function isdisjoint(U1::UnionSetArray, U2::UnionSetArray, witness::Bool=false)
    return _isdisjoint_union(U1, U2, witness)
end

"""
    isdisjoint(U::Universe, X::LazySet, [witness]::Bool=false)

Check whether a universe and another set do not intersect, and otherwise
optionally compute a witness.

### Input

- `U`       -- universe
- `X`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

`true` iff ``X ≠ ∅``.
"""
@commutative function isdisjoint(U::Universe, X::LazySet, witness::Bool=false)
    return _isdisjoint_universe(U, X, witness)
end

function _isdisjoint_universe(U::Universe, X::LazySet, witness)
    @assert dim(X) == dim(U) "the dimensions of the given sets should match, " *
                             "but they are $(dim(X)) and $(dim(U)), respectively"
    result = isempty(X)
    if result
        return _witness_result_empty(witness, true, U, X)
    else
        return witness ? (false, an_element(X)) : false
    end
end

# disambiguations
for ST in [:AbstractPolyhedron, :AbstractZonotope, :AbstractSingleton,
           :HalfSpace, :Hyperplane, :Line2D, :CartesianProductArray, :UnionSet,
           :UnionSetArray, :Complement]
    @eval @commutative function isdisjoint(U::Universe, X::($ST), witness::Bool=false)
        return _isdisjoint_universe(U, X, witness)
    end
end

function isdisjoint(U::Universe, ::Universe, witness::Bool=false)
    return witness ? (false, an_element(U)) : false
end

"""
    isdisjoint(C::Complement, X::LazySet, [witness]::Bool=false)

Check whether the complement of a set and another set do not intersect, and
otherwise optionally compute a witness.

### Input

- `C`       -- complement of a set
- `X`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

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
@commutative function isdisjoint(C::Complement, X::LazySet, witness::Bool=false)
    return _isdisjoint_complement(C, X, witness)
end

function _isdisjoint_complement(C, X, witness)
    return ⊆(X, C.X, witness)
end

# disambiguations
for ST in [:AbstractPolyhedron, :AbstractSingleton, :UnionSet, :UnionSetArray,
           :Hyperplane, :Line2D, :HalfSpace]
    @eval @commutative function isdisjoint(C::Complement, X::($ST), witness::Bool=false)
        return _isdisjoint_complement(C, X, witness)
    end
end

function isdisjoint(C1::Complement, C2::Complement, witness::Bool=false)
    return _isdisjoint_general(C1, C2, witness)
end

"""
    isdisjoint(cpa::CartesianProductArray, P::AbstractPolyhedron,
               [witness]::Bool=false)

Check whether a polytopic Cartesian product array and a polyhedral set do not
intersect, and otherwise optionally compute a witness.

### Input

- `cpa`     -- Cartesian products of a finite number of polytopes
- `P`       -- polyhedral set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``\\text{cpa} ∩ P = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``\\text{cpa} ∩ P = ∅``
  * `(false, v)` iff ``\\text{cpa} ∩ P ≠ ∅`` and ``v ∈ \\text{cpa} ∩ P``

### Algorithm

We first identify the blocks of `cpa` in which `P` is constrained.
Then we project `cpa` to those blocks and convert the result to an `HPolytope`
(or `HPolyhedron` if the set type is not known to be bounded) `Q`.
Finally we determine whether `Q` and the projected `P` intersect.
"""
@commutative function isdisjoint(cpa::CartesianProductArray,
                                 P::AbstractPolyhedron, witness::Bool=false)
    return _isdisjoint_cpa_polyhedron(cpa, P, witness)
end

function _isdisjoint_cpa_polyhedron(cpa::CartesianProductArray, P, witness)
    cpa_low_dim, vars, _block_structure = get_constrained_lowdimset(cpa, P)
    if !is_polyhedral(cpa_low_dim)
        error("a polyhedral set is required")
    end
    T = isconvextype(typeof(cpa_low_dim)) ? HPolytope : HPolyhedron
    hpoly_low_dim = T(constraints_list(cpa_low_dim))
    return isdisjoint(hpoly_low_dim, project(P, vars), witness)
end

# disambiguations
for ST in [:Hyperplane, :Line2D]
    @eval @commutative function isdisjoint(cpa::CartesianProductArray, X::($ST),
                                           witness::Bool=false)
        return _isdisjoint_cpa_polyhedron(cpa, X, witness)
    end
end

"""
    isdisjoint(X::CartesianProductArray, Y::CartesianProductArray,
               [witness]::Bool=false)

Check whether two Cartesian products of a finite number of sets with the same
block structure do not intersect, and otherwise optionally compute a witness.

### Input

- `X`       -- Cartesian products of a finite number of sets
- `Y`       -- Cartesian products of a finite number of sets
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ∩ Y = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``X ∩ Y = ∅``
  * `(false, v)` iff ``X ∩ Y ≠ ∅`` and ``v ∈ X ∩ Y``

### Notes

The implementation requires (and checks) that the Cartesian products have the
same block structure.

Witness production is currently not supported.
"""
function isdisjoint(X::CartesianProductArray, Y::CartesianProductArray,
                    witness::Bool=false)
    @assert same_block_structure(array(X), array(Y)) "block structure has to " *
                                                     "be identical"

    for i in eachindex(X.array)
        if isdisjoint(X.array[i], Y.array[i])
            return _witness_result_empty(witness, true, X, Y)
        end
    end
    if witness
        error("witness production is not supported yet")
    else
        return false
    end
end

"""
    isdisjoint(cpa::CartesianProductArray, H::AbstractHyperrectangle,
               [witness]::Bool=false)

Check whether a Cartesian product of a finite number of sets and a
hyperrectangular set do not intersect, and otherwise optionally compute a
witness.

### Input

- `cpa`     -- Cartesian product of a finite number of sets
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
@commutative function isdisjoint(cpa::CartesianProductArray,
                                 H::AbstractHyperrectangle,
                                 witness::Bool=false)
    N = promote_type(eltype(cpa), eltype(H))
    if witness
        w = zeros(N, dim(H))
    end
    n = dim(H)
    block_start = 1
    for bi in cpa
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

@commutative function isdisjoint(∅::EmptySet, X::LazySet, witness::Bool=false)
    return _isdisjoint_emptyset(∅, X, witness)
end

function _isdisjoint_emptyset(∅::EmptySet, X::LazySet, witness::Bool=false)
    @assert dim(∅) == dim(X) "the dimensions of the given sets should match, " *
                             "but they are $(dim(∅)) and $(dim(X)), respectively"
    return _witness_result_empty(witness, true, ∅, X)
end

# disambiguations
for ST in [:AbstractPolyhedron, :AbstractSingleton, :HalfSpace, :Hyperplane,
           :Line2D, :Universe, :Complement, :UnionSet, :UnionSetArray]
    @eval @commutative function isdisjoint(∅::EmptySet, X::($ST), witness::Bool=false)
        return _isdisjoint_emptyset(∅, X, witness)
    end
end

isdisjoint(∅1::EmptySet, ∅2::EmptySet, witness::Bool=false) = _isdisjoint_emptyset(∅1, ∅2, witness)

"""
    isdisjoint(L1::Line2D, L2::Line2D, [witness]::Bool=false)

Check whether two two-dimensional lines do not intersect, and otherwise
optionally compute a witness.

### Input

- `L1`      -- two-dimensional line
- `L2`      -- two-dimensional line
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``L1 ∩ L2 = ∅``
* If `witness` option is activated:
  * `(true, [])` iff ``L1 ∩ L2 = ∅``
  * `(false, v)` iff ``L1 ∩ L2 ≠ ∅`` and ``v ∈ L1 ∩ L2``
"""
function isdisjoint(L1::Line2D, L2::Line2D, witness::Bool=false)
    disjoint = _isdisjoint(L1, L2)
    if disjoint
        return _witness_result_empty(witness, true, L1, L2)
    end
    return witness ? (false, an_element(intersection(L1, L2))) : false
end

# the lines do not intersect <=> det is zero and they are not identical
function _isdisjoint(L1::Line2D, L2::Line2D)
    det = right_turn(L1.a, L2.a)
    disjoint = isapproxzero(det) && !isapprox(L1.b, L2.b)
    return disjoint
end

for ST in [:AbstractZonotope, :AbstractSingleton]
    @eval @commutative function isdisjoint(C::CartesianProduct{N,<:LazySet,<:Universe},
                                           Z::$(ST)) where {N}
        X = first(C)
        Zp = project(Z, 1:dim(X))
        return isdisjoint(X, Zp)
    end

    @eval @commutative function isdisjoint(C::CartesianProduct{N,<:Universe,<:LazySet},
                                           Z::$(ST)) where {N}
        Y = second(C)
        Zp = project(Z, (dim(first(C)) + 1):dim(C))
        return isdisjoint(Y, Zp)
    end

    # disambiguation
    @eval @commutative function isdisjoint(::CartesianProduct{N,<:Universe,<:Universe},
                                           Z::$(ST)) where {N}
        return false
    end
end

# see Proposition 8 in Wetzlinger, Kochdumper, Bak, Althoff: * Fully-automated
# verification of linear systems using inner- and outer-approximations of
# reachable sets*. 2022.
@commutative function isdisjoint(Z::AbstractZonotope, P::AbstractPolyhedron,
                                 witness::Bool=false; solver=nothing)
    n = dim(Z)
    @assert n == dim(P) "incompatible dimensions $(dim(Z)) and $(dim(P))"

    if n <= 2
        # this implementation is slower for low-dimensional sets
        return _isdisjoint_polyhedron(Z, P, witness; solver=solver)
    end

    N = promote_type(eltype(Z), eltype(P))
    c = center(Z)
    G = genmat(Z)
    C, d = tosimplehrep(P)
    p = size(G, 2)
    m = length(d)

    A = [C zeros(N, m, p);
         I(n) -G]
    b = vcat(d, c)
    obj = zeros(N, size(A, 2))

    lbounds = vcat(fill(-Inf, n), fill(-one(N), p))
    ubounds = vcat(fill(Inf, n), fill(one(N), p))
    sense = vcat(fill('<', m), fill('=', n))
    if isnothing(solver)
        solver = default_lp_solver(N)
    end

    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)

    if is_lp_optimal(lp.status)
        disjoint = false
    elseif is_lp_infeasible(lp.status)
        disjoint = true
    else
        return error("LP returned status $(lp.status) unexpectedly")
    end

    if disjoint
        return _witness_result_empty(witness, true, Z, P)
    elseif witness
        w = lp.sol[1:n]
        return (false, w)
    else
        return false
    end
end
