"""
    minkowski_sum(P::LazySet, Q::LazySet;
                  [backend]=nothing, [algorithm]=nothing, [prune]=true)

Compute the Minkowski sum of two polyhedral sets.

### Input

- `P`         -- set
- `Q`         -- set
- `backend`   -- (optional, default: `nothing`) polyhedral computations backend
- `algorithm` -- (optional, default: `nothing`) algorithm to eliminate
                 variables; available options are `Polyhedra.FourierMotzkin`,
                 `Polyhedra.BlockElimination`, and `Polyhedra.ProjectGenerators`
- `prune`     -- (optional, default: `true`) if `true`, apply a post-processing
                 to remove redundant constraints or vertices

### Output

In two dimensions, if the sets are polygons, the result is a `VPolygon`. In
higher dimensions, the result is an `HPolytope` if both `P` and `Q` are known to
be bounded by their types, and an `HPolyhedron` otherwise.

### Notes

This function requires that the list of constraints of both sets `P` and `Q` can
be obtained. After obtaining the respective lists of constraints, the
`minkowski_sum` method for polyhedral sets is used.
"""
function minkowski_sum(P::LazySet, Q::LazySet;
                       backend=nothing, algorithm=nothing, prune=true)
    n = dim(P)
    @assert n == dim(Q) "expected that the sets have the same dimension, " *
                        "but they are $n and $(dim(Q)) respectively"

    @assert is_polyhedral(P) && is_polyhedral(Q) "this function requires " *
                                                 "polyhedral sets; try overapproximating with an `HPolytope` or " *
                                                 "`HPolyhedron` first"

    if n == 2 && isboundedtype(typeof(P)) && isboundedtype(typeof(Q))
        # use vertex representation
        return _minkowski_sum_vpolygon(P, Q)
    end

    # use constraint representation
    res = _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
    if isboundedtype(typeof(P)) && isboundedtype(typeof(Q))
        return convert(HPolytope, res)
    else
        return res
    end
end

"""
    minkowski_sum(P::AbstractPolyhedron, Q::AbstractPolyhedron;
                  [backend]=nothing, [algorithm]=nothing, [prune]=true)

Compute the Minkowski sum of two polyhedra in constraint representation.

### Input

- `P`         -- polyhedron in constraint representation
- `Q`         -- polyhedron in constraint representation
- `backend`   -- (optional, default: `nothing`) polyhedral computations backend
- `algorithm` -- (optional, default: `nothing`) algorithm to eliminate
                 variables; available options are `Polyhedra.FourierMotzkin`,
                 `Polyhedra.BlockElimination`, and `Polyhedra.ProjectGenerators`
- `prune`     -- (optional, default: `true`) if `true`, apply a post-processing
                 to remove redundant constraints

### Output

A polyhedron in H-representation that corresponds to the Minkowski sum of `P`
and `Q`.

### Algorithm

This function implements the concrete Minkowski sum by projection and variable
elimination as detailed in [1]. The idea is that if we write ``P`` and ``Q`` in
*simple H-representation*, that is, ``P = \\{x ∈ ℝ^n : Ax ≤ b \\}``
and ``Q = \\{x ∈ ℝ^n : Cx ≤ d \\}``, then their Minkowski sum can be
seen as the projection onto the first ``n``-dimensional coordinates of the
polyhedron:
```math
    \\begin{pmatrix} 0 & A \\ C & -C \\end{pmatrix} \\binom{x}{y} ≤ \\binom{b}{d}
```
This is seen by noting that ``P ⊕ Q`` corresponds to the set of points
``x ∈ ℝ^n`` such that ``x = y + z`` with ``Ay ≤ b`` and ``Cz ≤ d``;
hence it follows that ``Ay ≤ b`` and ``C(x-y) ≤ d``, and the inequality
above follows by considering the ``2n``-dimensional space ``\\binom{x}{y}``.
The reduction from ``2n`` to ``n`` variables is performed using an elimination
algorithm as described next.

The elimination of variables depends on the polyhedra library `Polyhedra`, which
itself uses `CDDLib` for variable elimination. The available algorithms are:

- `Polyhedra.FourierMotzkin`    -- projection by computing the H-representation
                                   and applying the Fourier-Motzkin elimination
                                   algorithm to it

- `Polyhedra.BlockElimination`  -- projection by computing the H-representation
                                   and applying the block elimination algorithm
                                   to it

- `Polyhedra.ProjectGenerators` -- projection by computing the V-representation

[1] Kvasnica, Michal. "Minkowski addition of convex polytopes." (2005): 1-10.
"""
function minkowski_sum(P::AbstractPolyhedron, Q::AbstractPolyhedron;
                       backend=nothing, algorithm=nothing, prune=true)
    return _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
end

function minkowski_sum(P::AbstractPolytope, Q::AbstractPolytope;
                       backend=nothing, algorithm=nothing, prune=true)
    n = dim(P)
    @assert n == dim(Q) "expected that the sets have the same dimension, " *
                        "but they are $n and $(dim(Q)) respectively"

    if n == 2
        return _minkowski_sum_vpolygon(P, Q)
    end

    # fallback
    return _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
end

# common code before calling _minkowski_sum_hrep
function _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
    require(@__MODULE__, :Polyhedra; fun_name="minkowski_sum")
    require(@__MODULE__, :CDDLib; fun_name="minkowski_sum")

    A, b = tosimplehrep(P)
    C, d = tosimplehrep(Q)
    return _minkowski_sum_hrep(A, b, C, d; backend=backend, algorithm=algorithm,
                               prune=prune)
end

# This function computes the Minkowski sum of two polyhedra in simple
# H-representation, P = {x : Ax <= b} and Q = {x : Cx <= d}, using projection
# methods. See the documentation of `minkowski_sum` for details.
function _minkowski_sum_hrep(A::AbstractMatrix, b::AbstractVector,
                             C::AbstractMatrix, d::AbstractVector;
                             backend=nothing, algorithm=nothing, prune=true)
    if isnothing(backend)
        N = promote_type(eltype(A), eltype(b), eltype(C), eltype(d))
        backend = default_cddlib_backend(N)
    end

    if isnothing(algorithm)
        algorithm = Polyhedra.FourierMotzkin()
    elseif !(algorithm isa Polyhedra.EliminationAlgorithm)
        error("algorithm $algorithm is not a valid elimination algorithm; " *
              "choose among any of $(subtypes(Polyhedra.EliminationAlgorithm))")
    end

    mP, nP = size(A)
    mQ, nQ = size(C)
    E = [zeros(N, mP, nQ) A; C -C]
    f = [b; d]
    PQ = HPolyhedron(E, f)
    PQ_cdd = polyhedron(PQ; backend=backend)
    W_cdd = Polyhedra.eliminate(PQ_cdd, (nP + 1):(2nP), algorithm)
    W = convert(HPolyhedron, W_cdd)
    if prune
        success = remove_redundant_constraints!(W)
        if !success
            error("the constraints corresponding to the Minkowski sum of the " *
                  "given sets are infeasible")
        end
    end
    return W
end

"""
    minkowski_sum(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle)

Concrete Minkowski sum of a pair of hyperrectangular sets.

### Input

- `H1` -- hyperrectangular set
- `H2` -- hyperrectangular set

### Output

A `Hyperrectangle` corresponding to the Minkowski sum of `H1` and `H2`.

### Algorithm

The resulting hyperrectangle is obtained by summing up the centers and radii.
"""
function minkowski_sum(H1::AbstractHyperrectangle, H2::AbstractHyperrectangle)
    c = center(H1) + center(H2)
    r = radius_hyperrectangle(H1) + radius_hyperrectangle(H2)
    return Hyperrectangle(c, r)
end

"""
    minkowski_sum(Z1::AbstractZonotope, Z2::AbstractZonotope)

Concrete Minkowski sum of a pair of zonotopic sets.

### Input

- `Z1` -- zonotopic set
- `Z2` -- zonotopic set

### Output

A `Zonotope` corresponding to the Minkowski sum of `Z1` and `Z2`.

### Algorithm

The resulting zonotope is obtained by summing up the centers and concatenating
the generators of `Z1` and `Z2`.
"""
function minkowski_sum(Z1::AbstractZonotope, Z2::AbstractZonotope)
    cnew = center(Z1) + center(Z2)
    Gnew = hcat(genmat(Z1), genmat(Z2))
    return Zonotope(cnew, Gnew)
end

"""
    minkowski_sum(X::AbstractSingleton, Y::AbstractSingleton)

Concrete Minkowski sum of a pair of singletons.

### Input

- `X` -- singleton
- `Y` -- singleton

### Output

A singleton

### Algorithm

The singleton obtained by summing the elements in `X` and `Y`.
"""
function minkowski_sum(X::AbstractSingleton, Y::AbstractSingleton)
    @assert dim(X) == dim(Y) "expected that the singletons have the same " *
                             "dimension, but they are $(dim(X)) and $(dim(Y)) respectively"
    return Singleton(element(X) + element(Y))
end

function _minkowski_sum_vpolygon(P::LazySet, Q::LazySet)
    return minkowski_sum(convert(VPolygon, P), convert(VPolygon, Q))
end

function _minkowski_sum_vrep_2d(vlistP::AbstractVector{<:AbstractVector{N}},
                                vlistQ::AbstractVector{<:AbstractVector{N}}) where {N}
    mP = length(vlistP)
    mQ = length(vlistQ)
    if mP == 1 || mQ == 1
        return _minkowski_sum_vrep_2d_singleton(vlistP, vlistQ)
    end

    EAST = N[1, 0]
    ORIGIN = N[0, 0]
    R = fill(ORIGIN, mP + mQ)
    k = _σ_helper(EAST, vlistP)
    j = _σ_helper(EAST, vlistQ)

    i = 1
    while i <= size(R, 1)
        P₁, P₂ = vlistP[(k - 1) % mP + 1], vlistP[(k % mP + 1)]
        P₁P₂ = P₂ - P₁
        Q₁, Q₂ = vlistQ[(j - 1) % mQ + 1], vlistQ[(j % mQ + 1)]
        Q₁Q₂ = Q₂ - Q₁
        R[i] = P₁ + Q₁
        turn = right_turn(P₁P₂, Q₁Q₂, ORIGIN)
        if turn > 0
            k += 1
        elseif turn < 0
            j += 1
        else
            pop!(R)
            k += 1
            j += 1
        end
        i += 1
    end
    return R
end

# assume that at least one of the arguments has length 1
function _minkowski_sum_vrep_2d_singleton(vlistP::Vector{VT},
                                          vlistQ::Vector{VT}) where {N,VT<:AbstractVector{N}}
    mP = length(vlistP)
    mQ = length(vlistQ)

    if min(mP, mQ) != 1
        throw(ArgumentError("expected one argument to have only one vertex, " *
                            "but got $mP and $mQ respectively"))
    elseif mP != 1
        return _minkowski_sum_vrep_2d_singleton(vlistQ, vlistP)
    end

    # vlistP has only one vertex, mP == 1
    R = Vector{VT}(undef, mQ)
    @inbounds begin
        p = vlistP[1]
        for i in 1:mQ
            R[i] = p + vlistQ[i]
        end
    end
    return R
end

function _minkowski_sum_vrep_nd(vlist1::Vector{VT}, vlist2::Vector{VT};
                                apply_convex_hull::Bool=true, backend=nothing,
                                solver=nothing) where {N,VT<:AbstractVector{N}}
    n, m = length(vlist1), length(vlist2)
    Vout = Vector{VT}()
    sizehint!(Vout, n * m)
    for vi in vlist1
        for vj in vlist2
            push!(Vout, vi + vj)
        end
    end
    if apply_convex_hull
        if isnothing(backend)
            require(@__MODULE__, :Polyhedra; fun_name="minkowski_sum")
            backend = default_polyhedra_backend_nd(N)
            solver = default_lp_solver_polyhedra(N)
        end
        convex_hull!(Vout; backend=backend, solver=solver)
    end
    return Vout
end

"""
    minkowski_sum(PZ::DensePolynomialZonotope, Z::AbstractZonotope)

Compute the Minkowski sum of a polynomial zonotope and a zonotopic set.

### Input

- `PZ` -- polynomial zonotope
- `Z`  -- zonotopic set

## Output

A polynomial zonotope whose center is the sum of the centers of `PZ` and `Z`
and whose generators are the concatenation of the generators of `PZ` and `Z`.
"""
@commutative function minkowski_sum(PZ::DensePolynomialZonotope,
                                    Z::AbstractZonotope)
    c = PZ.c + center(Z)
    G = [PZ.G genmat(Z)]
    return DensePolynomialZonotope(c, PZ.E, PZ.F, G)
end

# ZeroSet is the neutral element (+ disambiguation)
for T in [:LazySet, :AbstractPolyhedron, :AbstractPolytope, :AbstractZonotope,
          :AbstractHyperrectangle, :AbstractSingleton, :DensePolynomialZonotope,
          :SparsePolynomialZonotope]
    @eval begin
        @commutative minkowski_sum(::ZeroSet, X::$T) = X
    end
end

# See Proposition 3.1.19 in [1].
# [1] Kochdumper. *Extensions of polynomial zonotopes and their application to
#     verification of cyber-physical systems.* PhD diss., TU Munich, 2022.
@commutative function minkowski_sum(PZ::SparsePolynomialZonotope, Z::AbstractZonotope)
    c = center(PZ) + center(Z)
    G = genmat_dep(PZ)
    GI = hcat(genmat_indep(PZ), genmat(Z))
    E = expmat(PZ)
    return SparsePolynomialZonotope(c, G, GI, E)
end

# Given two balls of the same p-norm, their Minkowski sum is again a p-norm ball,
# which can be seen as follows:
#
# σ_Bp(d) = c + r \\frac{d}{‖d‖_p}
#
#   ρ_Bp(d)
# = ⟨c + r \\frac{d}{‖d‖_p}, d⟩
# = ⟨c, d⟩ + \\frac{r}{‖d‖_p} * ⟨d, d⟩
#
#   ρ_Bp¹⊕Bp²(d)
# = ρ_Bp¹(d) + ρ_Bp²(d)
# = ⟨c¹, d⟩ + \\frac{r¹}{‖d‖_p} * ⟨d, d⟩ + ⟨c², d⟩ + \\frac{r²}{‖d‖_p} * ⟨d, d⟩
# = ⟨c¹ + c², d⟩ + \\frac{r¹ + r²}{‖d‖_p} * ⟨d, d⟩
# = ρ_Bp³(d)
#
# where Bp³ = Ball(center(Bp¹) + center(Bp²), radius(Bp¹) + radius(Bp²))
function minkowski_sum(B1::BT1, B2::BT2) where {BT<:AbstractBallp,BT1<:BT,BT2<:BT}
    p = ball_norm(B1)
    if ball_norm(B2) != p
        throw(ArgumentError("this method only applies to balls of the same norm"))
    end
    return Ballp(p, center(B1) + center(B2), radius_ball(B1) + radius_ball(B2))
end

for B in (:Ball1, :BallInf)
    @eval begin
        # see AbstractBallp method
        function minkowski_sum(B1::$B, B2::$B)
            return $B(center(B1) + center(B2), radius_ball(B1) + radius_ball(B2))
        end
    end
end
