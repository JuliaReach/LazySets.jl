export minkowski_sum

"""
    minkowski_sum(P::LazySet, Q::LazySet;
                  [backend]=nothing,
                  [algorithm]=nothing,
                  [prune]=true)

Concrete Minkowski sum for a pair of lazy sets using their constraint representation.

### Input

- `P`         -- lazy set
- `Q`         -- another lazy set
- `backend`   -- (optional, default: `nothing`) polyhedral computations backend
- `algorithm` -- (optional, default: `nothing`) algorithm to compute the elimination
                 of variables; available options are `Polyhedra.FourierMotzkin`,
                 `Polyhedra.BlockElimination`, and `Polyhedra.ProjectGenerators`
- `prune`     -- (optional, default: `true`) if `true`, apply a post-processing algorithm
                 to remove redundant constraints

### Output

In two dimensions the result is a `VPolygon`. In higher dimensions, the result
is an `HPolytope` if both `P` and `Q` are bounded, and an `HPolyhedron`
otherwise.

### Notes

This function requires that the list of constraints of both lazy sets `P` and
`Q` can be obtained. After obtaining the respective lists of constraints, the
`minkowski_sum` function for polyhedral sets is used. For details see
[`minkowski_sum(::VPolytope, ::VPolytope)`](@ref).

This method requires `Polyhedra` and `CDDLib`, so you have to do:

```julia
julia> using LazySets, Polyhedra, CDDLib

julia> ...

julia> minkowski_sum(P, Q)
```
"""
function minkowski_sum(P::LazySet, Q::LazySet;
                       backend=nothing,
                       algorithm=nothing,
                       prune=true)
    n = dim(P)
    @assert n == dim(Q) "expected that the sets have the same dimension, " *
                        "but they are $n and $(dim(Q)) respectively"

    if n == 2 && applicable(vertices_list, P) && applicable(vertices_list, Q) &&
                 isboundedtype(typeof(P)) && isboundedtype(typeof(Q))
        Pv = vertices_list(P)
        if length(Pv) > 1
            Pv = _convex_hull_2d_preprocess!(copy(Pv))
        end
        Qv = vertices_list(Q)
        if length(Qv) > 1
            Qv = _convex_hull_2d_preprocess!(copy(Qv))
        end
        R = _minkowski_sum_vrep_2d(Pv, Qv)
        return VPolygon(R)
        # return _minkowski_sum_vpolygon(P, Q) # crashes, see JuliaLang#41561
    end

    @assert applicable(constraints_list, P) &&
        applicable(constraints_list, Q) "this function requires that the " *
        "list of constraints is available for both arguments; try " *
        "overapproximating with an `HPolytope` first"

    res = _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
    if isbounded(P) && isbounded(Q)
        return convert(HPolytope, res)
    else
        return res
    end
end

"""
    minkowski_sum(P::AbstractPolyhedron, Q::AbstractPolyhedron;
                  [backend]=nothing,
                  [algorithm]=nothing,
                  [prune]=true)

Compute the Minkowski sum between two polyhedra in constraint representation.

### Input

- `P`         -- polyhedron in constraint representation
- `Q`         -- another polyhedron in constraint representation
- `backend`   -- (optional, default: `nothing`) polyhedral computations backend
- `algorithm` -- (optional, default: `nothing`) algorithm to compute the elimination
                 of variables; available options are `Polyhedra.FourierMotzkin`,
                 `Polyhedra.BlockElimination`, and `Polyhedra.ProjectGenerators`
- `prune`     -- (optional, default: `true`) if `true`, apply a post-processing algorithm
                 to remove redundant constraints

### Output

A polyhedron in H-representation that corresponds to the Minkowski sum of `P` and `Q`.

### Notes

This method requires `Polyhedra` and `CDDLib`, so you have to do:

```julia
julia> using LazySets, Polyhedra, CDDLib

julia> ...

julia> minkowski_sum(P, Q)
```

### Algorithm

This function implements the concrete Minkowski sum by projection and variable
elimination as detailed in [1]. The idea is that if we write ``P`` and ``Q`` in
*simple H-representation*, that is, ``P = \\{x ∈ \\mathbb{R}^n : Ax ≤ b \\}``
and ``Q = \\{x ∈ \\mathbb{R}^n : Cx ≤ d \\}``, then their Minkowski sum can be
seen as the projection onto the first ``n``-dimensional coordinates of the polyhedron

```math
    \\begin{pmatrix} 0 & A \\ C & -C \\end{pmatrix} \\binom{x}{y} ≤ \binom{b}{d}
```
This is seen by noting that ``P ⊕ Q`` corresponds to the set of points
``x ∈ \\mathbb{R}^n`` such that ``x = y + z`` with ``Ay ≤ b`` and ``Cz ≤ d``;
hence it follows that ``Ay ≤ b`` and ``C(x-y) ≤ d``, and the inequality displayed
above follows by considering the ``2n``-dimensional space ``\\binom{x}{y}``.
The reduction from ``2n`` to ``n`` variables is performed using an elimination
algorithm as described next.

The elimination of variables depends on the concrete polyhedra library `Polyhedra`,
which itself uses `CDDLib` for variable elimination. The available algorithms are:

- `Polyhedra.FourierMotzkin`   -- computation of the projection by computing the
                                  H-representation and applying the Fourier-Motzkin
                                  elimination algorithm to it

- `Polyhedra.BlockElimination` -- computation of the projection by computing the
                                  H-representation and applying the block elimination
                                  algorithm to it

- `Polyhedra.ProjectGenerators` -- computation of the projection by computing the
                                   V-representation

[1] Kvasnica, Michal. "Minkowski addition of convex polytopes." (2005): 1-10.
"""
function minkowski_sum(P::AbstractPolyhedron, Q::AbstractPolyhedron;
                       backend=nothing,
                       algorithm=nothing,
                       prune=true)
    return _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
end

function minkowski_sum(P::AbstractPolytope, Q::AbstractPolytope;
                       backend=nothing,
                       algorithm=nothing,
                       prune=true)
    n = dim(P)
    @assert n == dim(Q) "expected that the sets have the same dimension, " *
                        "but they are $n and $(dim(Q)) respectively"

    if n == 2
        return _minkowski_sum_vpolygon(P, Q)
    end

    # fallback
    return _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
end

function minkowski_sum(P::HPolytope, Q::HPolytope;
                       backend=nothing,
                       algorithm=nothing,
                       prune=true)
    res = _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
    return convert(HPolytope, res)
end

# common code before calling _minkowski_sum_hrep
function _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
    require(:Polyhedra; fun_name="minkowski_sum")
    require(:CDDLib; fun_name="minkowski_sum")

    A, b = tosimplehrep(P)
    C, d = tosimplehrep(Q)
    return _minkowski_sum_hrep(A, b, C, d, backend=backend, algorithm=algorithm,
                               prune=prune)
end

# This function computes the concrete Minkowski sum between two polyhedra in
# simple H-representation,
# P = {x : Ax <= b} and Q = {x : Cx <= d}
# using the projection methods. See the documentation of `minkowski_sum` for details.
function _minkowski_sum_hrep(A::AbstractMatrix, b::AbstractVector,
                             C::AbstractMatrix, d::AbstractVector;
                             backend=nothing,
                             algorithm=nothing,
                             prune=true)

    if backend == nothing
        N = promote_type(eltype(A), eltype(b), eltype(C), eltype(d))
        backend = default_cddlib_backend(N)
    end

    if algorithm == nothing
        algorithm = Polyhedra.FourierMotzkin()
    elseif !(algorithm <: EliminationAlgorithm)
        error("the algorithm $algorithm is not a valid elimination algorithm;
              choose among any of $(subtypes(Polyhedra.EliminationAlgorithm))")
    end

    mP, nP = size(A)
    mQ, nQ = size(C)
    E = [zeros(N, mP, nQ) A; C -C]
    f = [b; d]
    PQ = HPolyhedron(E, f)
    PQ_cdd = polyhedron(PQ, backend=backend)
    W_cdd = Polyhedra.eliminate(PQ_cdd, nP+1:2nP, algorithm)
    W = convert(HPolyhedron, W_cdd)
    if prune
        success = remove_redundant_constraints!(W)
        if !success
            error("the constraints corresponding to the minkowski sum of the given " *
                  "sets are infeasible")
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

A `Hyperrectangle` corresponding to the concrete Minkowski sum of `H1` and `H2`.

### Algorithm

The resulting hyperrectangle is obtained by summing up the centers and
radiuses of `H1` and `H2`.
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

A `Zonotope` corresponding to the concrete Minkowski sum of `Z1` and `Z2`.

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
    @assert dim(X) == dim(Y) "expected that the singletons have the same dimension, " *
                "but they are $(dim(X)) and $(dim(Y)) respectively"
    return Singleton(element(X) + element(Y))
end

"""
    minkowski_sum(x::Interval, y::Interval)

Concrete Minkowski sum of a pair of intervals.

### Input

- `x` -- hyperrectangular set
- `y` -- hyperrectangular set

### Output

An `Interval` corresponding to the concrete Minkowski sum of `x` and `y`.

### Algorithm

The function takes the sum of `x` and `y` following the rules of interval
arithmetic.
"""
function minkowski_sum(x::Interval, y::Interval)
    return Interval(x.dat + y.dat)
end

"""
    minkowski_sum(P::VPolygon, Q::VPolygon)

The Minkowski Sum of two polygon in vertex representation.

### Input

- `P` -- polygon in vertex representation
- `Q` -- another polygon in vertex representation

### Output

A polygon in vertex representation.

### Algorithm

We treat each edge of the polygons as a vector, attaching them in polar order
(attaching the tail of the next vector to the head of the previous vector). The
resulting polygonal chain will be a polygon, which is the Minkowski sum of the
given polygons. This algorithm assumes that the vertices of P and Q are sorted
in counter-clockwise fashion and has linear complexity O(m+n) where m and n are
the number of vertices of P and Q respectively.
"""
function minkowski_sum(P::VPolygon, Q::VPolygon)
    R = _minkowski_sum_vrep_2d(P.vertices, Q.vertices)
    return VPolygon(R)
end

function _minkowski_sum_vpolygon(P, Q)
    return minkowski_sum(convert(VPolygon, P), convert(VPolygon, Q))
end

function _minkowski_sum_vrep_2d(vlistP::Vector{VT},
                                vlistQ::Vector{VT}) where {N, VT<:AbstractVector{N}}
    mP = length(vlistP)
    mQ = length(vlistQ)
    if mP == 1 || mQ == 1
        return _minkowski_sum_vrep_2d_singleton(vlistP, vlistQ)
    end

    EAST = N[1, 0]
    ORIGIN = N[0, 0]
    k = _σ_helper(EAST, vlistP)
    j = _σ_helper(EAST, vlistQ)
    R = Vector{VT}(undef, mP+mQ)
    fill!(R, ORIGIN)

    i = 1
    while i <= size(R, 1)
        P₁, P₂ = vlistP[(k-1)%mP + 1], vlistP[(k%mP + 1)]
        P₁P₂ = P₂ - P₁
        Q₁, Q₂ = vlistQ[(j-1)%mQ + 1], vlistQ[(j%mQ + 1)]
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
                                          vlistQ::Vector{VT}) where {N, VT<:AbstractVector{N}}
    mP = length(vlistP)
    mQ = length(vlistQ)

    if min(mP, mQ) != 1
        throw(ArgumentError("expected one argument to have only one vertex, got $mP and $mQ respectively"))
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

"""
    minkowski_sum(P1::VPolytope, P2::VPolytope;
                  [apply_convex_hull]=true,
                  [backend]=nothing,
                  [solver]=nothing)

Compute the Minkowski sum between two polytopes in vertex representation.

### Input

- `P1`                -- polytope
- `P2`                -- another polytope
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         pairwise sums using a convex hull algorithm
- `backend`           -- (optional, default: `nothing`) the backend for
                         polyhedral computations used to post-process with a
                         convex hull; see `default_polyhedra_backend(P1)`
- `solver`            -- (optional, default: `nothing`) the backend used to
                         solve the linear program; see
                         `default_lp_solver_polyhedra(N)`

### Output

A new polytope in vertex representation whose vertices are the convex hull of
the sum of all possible sums of vertices of `P1` and `P2`.
"""
function minkowski_sum(P1::VPolytope, P2::VPolytope;
                       apply_convex_hull::Bool=true,
                       backend=nothing,
                       solver=nothing)

    @assert dim(P1) == dim(P2) "cannot compute the Minkowski sum between a polyotope " *
        "of dimension $(dim(P1)) and a polytope of dimension $((dim(P2)))"

    vlist1 = _vertices_list(P1, backend)
    vlist2 = _vertices_list(P2, backend)
    Vout = _minkowski_sum_vrep_nd(vlist1, vlist2, apply_convex_hull=apply_convex_hull, backend=backend, solver=solver)
    return VPolytope(Vout)
end

function _minkowski_sum_vrep_nd(vlist1::Vector{VT}, vlist2::Vector{VT};
                                apply_convex_hull::Bool=true,
                                backend=nothing,
                                solver=nothing) where {N, VT<:AbstractVector{N}}
    n, m = length(vlist1), length(vlist2)
    Vout = Vector{VT}()
    sizehint!(Vout, n * m)
    for vi in vlist1
        for vj in vlist2
            push!(Vout, vi + vj)
        end
    end
    if apply_convex_hull
        if backend == nothing
            require(:Polyhedra; fun_name="minkowski_sum")
            backend = default_polyhedra_backend_nd(N)
            solver = default_lp_solver_polyhedra(N)
        end
        convex_hull!(Vout, backend=backend, solver=solver)
    end
    return Vout
end

"""
    minkowski_sum(PZ::PolynomialZonotope, Z::AbstractZonotope)

Return the Minkowski sum of a polynomial zonotope and a usual zonotopic set.

### Input

- `PZ` -- polynomial zonotope
- `Z`  -- usual zonotopic set

## Output

A polynomial zonotope whose center is the sum of the centers of `PZ` and `Z`
and whose generators are the concatenation of the generators of `PZ` and `Z`.
"""
function minkowski_sum(PZ::PolynomialZonotope, Z::AbstractZonotope)
    c = PZ.c + center(Z)
    G = [PZ.G genmat(Z)]
    return PolynomialZonotope(c, PZ.E, PZ.F, G)
end

# symmetric method
minkowski_sum(Z::AbstractZonotope, PZ::PolynomialZonotope) = minkowski_sum(PZ, Z)

minkowski_sum(X::LazySet, ::ZeroSet) = X
minkowski_sum(::ZeroSet, X::LazySet) = X
minkowski_sum(Z::ZeroSet, ::ZeroSet) = Z

# disambiguation
minkowski_sum(::ZeroSet, P::AbstractPolyhedron) = P
minkowski_sum(P::AbstractPolyhedron, ::ZeroSet) = P
minkowski_sum(P::AbstractPolytope, ::ZeroSet) = P
minkowski_sum(::ZeroSet, P::AbstractPolytope) = P
minkowski_sum(::ZeroSet, Z::AbstractZonotope) = Z
minkowski_sum(Z::AbstractZonotope, ::ZeroSet) = Z
minkowski_sum(::ZeroSet, H::AbstractHyperrectangle) = H
minkowski_sum(H::AbstractHyperrectangle, ::ZeroSet) = H
minkowski_sum(X::AbstractSingleton, Y::ZeroSet) = X
minkowski_sum(::ZeroSet, X::AbstractSingleton) = X
