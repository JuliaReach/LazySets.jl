"""
# Extended help

    isdisjoint(X::LazySet, Y::LazySet, [witness]::Bool=false)

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
# Extended help

    isdisjoint(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle,
               [witness]::Bool=false)

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

function isdisjoint(S1::AbstractSingleton, S2::AbstractSingleton,
                    witness::Bool=false)
    s1 = element(S1)
    empty_intersection = !isapprox(s1, element(S2))
    return _witness_result_empty(witness, empty_intersection, S1, S2, s1)
end

# common code for singletons
function _isdisjoint_singleton(S::AbstractSingleton, X::LazySet,
                               witness::Bool=false)
    s = element(S)
    empty_intersection = s ∉ X
    return _witness_result_empty(witness, empty_intersection, S, X, s)
end

"""
# Extended help

    isdisjoint(X::LazySet, S::AbstractSingleton, [witness]::Bool=false)

### Algorithm

Let ``S = \\{s\\}``. Then ``S ∩ X = ∅`` iff ``s ∉ X``.
"""
@commutative function isdisjoint(X::LazySet, S::AbstractSingleton,
                                 witness::Bool=false)
    return _isdisjoint_singleton(S, X, witness)
end

"""
# Extended help

    isdisjoint(Z::AbstractZonotope, H::Hyperplane, [witness]::Bool=false)

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
# Extended help

    isdisjoint(Z1::AbstractZonotope, Z2::AbstractZonotope,
               [witness]::Bool=false; [solver]=nothing)

### Input

- `solver`  -- (optional, default: `nothing`) the backend used to solve the
               linear program

### Algorithm

The algorithm is taken from [GuibasNZ03](@citet).

``Z1 ∩ Z2 = ∅`` iff ``c_1 - c_2 ∉ Z(0, (g_1, g_2))`` where ``c_i`` and ``g_i``
are the center and generators of zonotope `Zi` and ``Z(c, g)`` represents the
zonotope with center ``c`` and generators ``g``.
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
# Extended help

    isdisjoint(X::LazySet, hp::Hyperplane, [witness]::Bool=false)

### Notes

This implementation assumes that the set `X` is convex.

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
@commutative function isdisjoint(X::LazySet, hp::Hyperplane, witness::Bool=false)
    return _isdisjoint_hyperplane(hp, X, witness)
end

@commutative function isdisjoint(X::LazySet, L::Line2D, witness::Bool=false)
    return _isdisjoint_hyperplane(L, X, witness)
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
# Extended help

    isdisjoint(X::LazySet, hs::HalfSpace, [witness]::Bool=false)

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
# Extended help

    isdisjoint(P::AbstractPolyhedron, X::LazySet, [witness]::Bool=false;
               [solver]=nothing, [algorithm]="exact")

### Input

- `solver`    -- (optional, default: `nothing`) the backend used to solve the
                 linear program
- `algorithm` -- (optional, default: `"exact"`) algorithm keyword, one of:
                 * `"exact" (exact, uses a feasibility LP)
                 * `"sufficient" (sufficient, uses half-space checks)

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
        if !ispolyhedral(X)
            error("this algorithm requires a polyhedral input")
        end
        clist_P = _normal_Vector(P) # TODO
        clist_X = _normal_Vector(X) # TODO
        if isnothing(solver)
            solver = default_lp_solver(N)
        end
        return _isinfeasible([clist_P; clist_X], witness; solver=solver)
    else
        error("algorithm $algorithm unknown")
    end
end

function isdisjoint(U1::UnionSet, U2::UnionSet, witness::Bool=false)
    return _isdisjoint_union(U1, U2, witness)
end

@commutative function isdisjoint(U::UnionSet, X::LazySet, witness::Bool=false)
    return _isdisjoint_union(U, X, witness)
end

function isdisjoint(U1::UnionSetArray, U2::UnionSetArray, witness::Bool=false)
    return _isdisjoint_union(U1, U2, witness)
end

@commutative function isdisjoint(U::UnionSetArray, X::LazySet, witness::Bool=false)
    return _isdisjoint_union(U, X, witness)
end

function _isdisjoint_union(cup::Union{UnionSet,UnionSetArray}, X::LazySet, witness::Bool=false)
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
    N = promote_type(eltype(cup), eltype(X))
    return _witness_result_empty(witness, true, N)
end

@commutative function isdisjoint(U::Universe, X::LazySet, witness::Bool=false)
    return _isdisjoint_universe(U, X, witness)
end

function isdisjoint(C1::Complement, C2::Complement, witness::Bool=false)
    return _isdisjoint_general(C1, C2, witness)
end

"""
# Extended help

    isdisjoint(C::Complement, X::LazySet, [witness]::Bool=false)

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

"""
# Extended help

    isdisjoint(cpa::CartesianProductArray, P::AbstractPolyhedron,
               [witness]::Bool=false)

### Notes

This implementation assumes that the sets in the Cartesian product `cpa` are polyhedral.

- `cpa`     -- Cartesian products of a finite number of polytopes

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
    if !ispolyhedral(cpa_low_dim)
        error("a polyhedral set is required")
    end
    T = isconvextype(typeof(cpa_low_dim)) ? HPolytope : HPolyhedron
    hpoly_low_dim = T(constraints_list(cpa_low_dim))
    return isdisjoint(hpoly_low_dim, project(P, vars), witness)
end

"""
# Extended help

    isdisjoint(X::CartesianProductArray, Y::CartesianProductArray,
               [witness]::Bool=false)

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
# Extended help

    isdisjoint(cpa::CartesianProductArray, H::AbstractHyperrectangle,
               [witness]::Bool=false)

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

@validate_commutative function isdisjoint(∅::EmptySet, X::LazySet, witness::Bool=false)
    return _isdisjoint_emptyset(∅, X, witness)
end

for T in [:AbstractZonotope, :AbstractSingleton]
    @eval @commutative function isdisjoint(C::CartesianProduct{N,<:LazySet,<:Universe},
                                           Z::$(T)) where {N}
        X = first(C)
        Zp = project(Z, 1:dim(X))
        return isdisjoint(X, Zp)
    end

    @eval @commutative function isdisjoint(C::CartesianProduct{N,<:Universe,<:LazySet},
                                           Z::$(T)) where {N}
        Y = second(C)
        Zp = project(Z, (dim(first(C)) + 1):dim(C))
        return isdisjoint(Y, Zp)
    end

    # disambiguation
    @eval @commutative function isdisjoint(::CartesianProduct{N,<:Universe,<:Universe},
                                           Z::$(T)) where {N}
        return false
    end
end

# See [WetzlingerKBA23; Proposition 8](@citet).
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

# ============== #
# disambiguation #
# ============== #

for T in (:AbstractPolyhedron, :AbstractZonotope, :AbstractHyperrectangle,
          :Hyperplane, :Line2D, :HalfSpace, :CartesianProductArray, :UnionSet,
          :UnionSetArray)
    @eval @commutative function isdisjoint(X::($T), S::AbstractSingleton, witness::Bool=false)
        return _isdisjoint_singleton(S, X, witness)
    end
end

for T in (:Hyperplane, :Line2D)
    @eval @commutative function isdisjoint(X::AbstractPolyhedron, Y::$(T), witness::Bool=false)
        return _isdisjoint_hyperplane(Y, X, witness)
    end
end

@commutative function isdisjoint(hp::Hyperplane, L::Line2D, witness::Bool=false)
    return _isdisjoint_hyperplane_hyperplane(hp, L, witness)
end

for T in (:AbstractPolyhedron, :AbstractZonotope, :Hyperplane, :Line2D, :CartesianProductArray)
    @eval @commutative function isdisjoint(X::($T), H::HalfSpace, witness::Bool=false)
        return _isdisjoint_halfspace(H, X, witness)
    end
end

function isdisjoint(X::AbstractPolyhedron, P::AbstractPolyhedron, witness::Bool=false;
                    solver=nothing, algorithm="exact")
    return _isdisjoint_polyhedron(P, X, witness)
end

for TU in (:UnionSet, :UnionSetArray)
    for T in (:AbstractPolyhedron, :Hyperplane, :Line2D, :HalfSpace)
        @eval @commutative function isdisjoint(U::($TU), X::($T), witness::Bool=false)
            return _isdisjoint_union(U, X, witness)
        end
    end
end

@commutative function isdisjoint(U1::UnionSet, U2::UnionSetArray, witness::Bool=false)
    return _isdisjoint_union(U1, U2, witness)
end

for T in (:AbstractPolyhedron, :AbstractZonotope, :AbstractSingleton,
          :HalfSpace, :Hyperplane, :Line2D, :CartesianProductArray, :UnionSet,
          :UnionSetArray, :Complement)
    @eval @commutative function isdisjoint(U::Universe, X::($T), witness::Bool=false)
        return _isdisjoint_universe(U, X, witness)
    end
end

for T in (:AbstractPolyhedron, :AbstractSingleton, :UnionSet, :UnionSetArray,
          :Hyperplane, :Line2D, :HalfSpace)
    @eval @commutative function isdisjoint(C::Complement, X::($T), witness::Bool=false)
        return _isdisjoint_complement(C, X, witness)
    end
end

for T in (:Hyperplane, :Line2D)
    @eval @commutative function isdisjoint(cpa::CartesianProductArray, X::($T), witness::Bool=false)
        return _isdisjoint_cpa_polyhedron(cpa, X, witness)
    end
end

for T in (:AbstractPolyhedron, :AbstractSingleton, :HalfSpace, :Hyperplane,
          :Line2D, :Universe, :Complement, :UnionSet, :UnionSetArray)
    @eval @validate_commutative function isdisjoint(∅::EmptySet, X::($T), witness::Bool=false)
        return _isdisjoint_emptyset(∅, X, witness)
    end
end
