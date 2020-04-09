export minkowski_sum

"""
    minkowski_sum(P::LazySet{N}, Q::LazySet{N};
                  [backend]=nothing,
                  [algorithm]=nothing,
                  [prune]=true) where {N<:Real}

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

An `HPolytope` that corresponds to the Minkowski sum of `P` and `Q` if both `P`
and `Q` are bounded; otherwise an `HPolyhedron`.

### Notes

This function requires that the list of constraints of both lazy sets `P` and
`Q` can be obtained. After obtaining the respective lists of constraints, the
`minkowski_sum` fucntion for polyhedral sets is used. For details see
[`minkowski_sum(::VPolytope, ::VPolytope)`](@ref).

This method requires `Polyhedra` and `CDDLib`, so you have to do:

```julia
julia> using LazySets, Polyhedra, CDDLib

julia> ...

julia> minkowski_sum(P, Q)
```
"""
function minkowski_sum(P::LazySet{N}, Q::LazySet{N};
                       backend=nothing,
                       algorithm=nothing,
                       prune=true) where {N<:Real}
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
    minkowski_sum(P::AbstractPolyhedron{N}, Q::AbstractPolyhedron{N};
                  [backend]=nothing,
                  [algorithm]=nothing,
                  [prune]=true) where {N<:Real}

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
function minkowski_sum(P::AbstractPolyhedron{N}, Q::AbstractPolyhedron{N};
                       backend=nothing,
                       algorithm=nothing,
                       prune=true) where {N<:Real}
    return _minkowski_sum_hrep_preprocess(P, Q, backend, algorithm, prune)
end

function minkowski_sum(P::HPolytope{N}, Q::HPolytope{N};
                       backend=nothing,
                       algorithm=nothing,
                       prune=true) where {N<:Real}
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
function _minkowski_sum_hrep(A::AbstractMatrix{N}, b::AbstractVector{N},
                             C::AbstractMatrix{N}, d::AbstractVector{N};
                             backend=nothing,
                             algorithm=nothing,
                             prune=true) where {N<:Real}

    if backend == nothing
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
    minkowski_sum(H1::AbstractHyperrectangle{N}, H2::AbstractHyperrectangle{N})
        where {N<:Real}

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
function minkowski_sum(H1::AbstractHyperrectangle{N},
                       H2::AbstractHyperrectangle{N}) where {N<:Real}
    return Hyperrectangle(center(H1) + center(H2),
                          radius_hyperrectangle(H1) + radius_hyperrectangle(H2))
end

"""
    minkowski_sum(Z1::AbstractZonotope{N}, Z2::AbstractZonotope{N}) where {N<:Real}

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
function minkowski_sum(Z1::AbstractZonotope{N}, Z2::AbstractZonotope{N}) where {N<:Real}
    cnew = center(Z1) + center(Z2)
    Gnew = hcat(genmat(Z1), genmat(Z2))
    return Zonotope(cnew, Gnew)
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
