export constrained_dimensions,
       remove_redundant_constraints,
       remove_redundant_constraints!,
       addconstraint!,
       ishyperplanar

ispolyhedraltype(::Type{<:AbstractPolyhedron}) = true

"""
# Extended help

    in(x::AbstractVector, P::AbstractPolyhedron)

### Algorithm

This implementation checks if the point lies inside each defining half-space.
"""
@validate function in(x::AbstractVector, P::AbstractPolyhedron)
    for c in constraints_list(P)
        if !_leq(dot(c.a, x), c.b)
            return false
        end
    end
    return true
end

"""
# Extended help

    isuniversal(P::AbstractPolyhedron, [witness]::Bool=false)

### Algorithm

`P` is universal iff it has no constraints.

A witness is produced using `isuniversal(H)` where `H` is the first linear
constraint of `P`.
"""
function isuniversal(P::AbstractPolyhedron, witness::Bool=false)
    c = constraints(P)
    if isempty(c)
        return _witness_result_empty(witness, true, eltype(P))
    else
        return witness ? isuniversal(first(c), true) : false
    end
end

struct LinearMapInverse{T,MT<:AbstractMatrix{T}} <: AbstractLinearMapAlgorithm
    inverse::MT
end

struct LinearMapInverseRight <: AbstractLinearMapAlgorithm end

struct LinearMapLift <: AbstractLinearMapAlgorithm end

struct LinearMapElimination{T,S} <: AbstractLinearMapAlgorithm
    backend::T
    method::S
end

struct LinearMapVRep{T} <: AbstractLinearMapAlgorithm
    backend::T
end

function _check_algorithm_applies(M::AbstractMatrix, P::LazySet,
                                  ::Type{LinearMapInverse};
                                  cond_tol=DEFAULT_COND_TOL,
                                  throw_error=false)
    inv_condition = issquare(M) && isinvertible(M; cond_tol=cond_tol)
    if !inv_condition
        throw_error && throw(ArgumentError("algorithm \"inverse\" requires " *
                                           "an invertible matrix"))
        return false
    end

    dense_condition = !issparse(M)
    if !dense_condition
        throw_error && throw(ArgumentError("the inverse of a sparse matrix " *
                                           "is not available; either convert your matrix to a dense matrix " *
                                           "with `Matrix(M)` or try the \"inverse_right\" algorithm"))
        return false
    end
    return true
end

function _check_algorithm_applies(M::AbstractMatrix, P::LazySet,
                                  ::Type{LinearMapInverseRight};
                                  cond_tol=DEFAULT_COND_TOL,
                                  throw_error=false)
    inv_condition = issquare(M) && isinvertible(M; cond_tol=cond_tol)
    if !inv_condition
        throw_error && throw(ArgumentError("algorithm \"inverse_right\" " *
                                           "requires an invertible matrix"))
        return false
    end
    return true
end

function _check_algorithm_applies(M::AbstractMatrix, P::LazySet,
                                  ::Type{LinearMapLift};
                                  throw_error=false)
    m, n = size(M)
    size_condition = m > n
    if !size_condition
        throw_error && throw(ArgumentError("algorithm \"lift\" requires that " *
                                           "the number of rows of the linear map is greater than the number " *
                                           "of columns, but they are of size $m and $n respectively"))
        return false
    end

    # rank condition
    r = rank(M)
    if r != n
        throw_error && throw(ArgumentError("the rank of the given matrix is " *
                                           "$r, but the algorithm \"lift\" requires it to be $n"))
        return false
    end

    return true
end

function _check_algorithm_applies(M::AbstractMatrix, P::LazySet,
                                  ::Type{LinearMapVRep};
                                  throw_error=false)
    if !isboundedtype(typeof(P))
        throw_error && throw(ArgumentError("algorithm \"vrep\" requires a " *
                                           "polytope and its list of vertices, but the set is a $(typeof(P))"))
        return false
    end
    return true
end

function _get_elimination_instance(N, backend, elimination_method)
    require(@__MODULE__, :Polyhedra; fun_name="linear_map with elimination")
    if isnothing(backend)
        require(@__MODULE__, :CDDLib; fun_name="linear_map with elimination")
        backend = default_cddlib_backend(N)
    end
    if isnothing(elimination_method)
        elimination_method = Polyhedra.BlockElimination()
    end
    return LinearMapElimination(backend, elimination_method)
end

function _default_linear_map_algorithm(M::AbstractMatrix, P::LazySet;
                                       cond_tol=DEFAULT_COND_TOL,
                                       backend=nothing,
                                       elimination_method=nothing)
    if _check_algorithm_applies(M, P, LinearMapInverse; cond_tol=cond_tol)
        algo = LinearMapInverse(inv(M))
    elseif _check_algorithm_applies(M, P, LinearMapLift)
        algo = LinearMapLift()
    else
        N = promote_type(eltype(M), eltype(P))
        algo = _get_elimination_instance(N, backend, elimination_method)
    end
    return algo
end

"""
    _linear_map_polyhedron(M::AbstractMatrix,
                           P::LazySet;
                           [algorithm]::Union{String, Nothing}=nothing,
                           [check_invertibility]::Bool=true,
                           [cond_tol]::Number=DEFAULT_COND_TOL,
                           [inverse]::Union{AbstractMatrix{N}, Nothing}=nothing,
                           [backend]=nothing,
                           [elimination_method]=nothing)

Concrete linear map of a polyhedral set.

### Input

- `M`         -- matrix
- `P`         -- polyhedral set
- `algorithm` -- (optional; default: `nothing`) algorithm to be used; for the
                 description see the Algorithm section below; possible choices
                 are:

    - `"inverse"`, alias: `"inv"`
    - `"inverse_right"`, alias: `"inv_right"`
    - `"elimination"`, alias: `"elim"`
    - `"lift"`
    - `"vrep"`
    - `"vrep_chull"`

- `check_invertibility` -- (optional, default: `true`) if `true` check whether
                           the given matrix `M` is invertible; set to `false`
                           only if you know that `M` is invertible
- `cond_tol`  -- (optional; default: `DEFAULT_COND_TOL`) tolerance of matrix
                 condition (used to check whether the matrix is invertible)
- `inverse`   -- (optional; default: `nothing`) matrix inverse `M⁻¹`; use this
                 option if you have already computed the inverse matrix of `M`
- `backend`   -- (optional: default: `nothing`) polyhedra backend
- `elimination_method`  -- (optional: default: `nothing`) elimination method for
                           the `"elimination"` algorithm

### Output

The type of the result is "as close as possible" to the the type of `P`.
Let `(m, n)` be the size of `M`, where `m ≠ n` is allowed for rectangular maps.

To fix the type of the output to something different than the default value,
consider post-processing the result of this function with a call to a suitable
`convert` method.

In particular, the output depends on the type of `P`, on `m`, and the algorithm
that was used:

- If the vertex-based approach was used:

    - If `P` is a `VPolygon` and `m = 2` then the output is a `VPolygon`.
    - If `P` is a `VPolytope` then the output is a `VPolytope`.
    - Otherwise the output is an `Interval` if `m = 1`, a `VPolygon` if `m = 2`,
      and a `VPolytope` in other cases.

- If the invertibility criterion was used:

    - The types of `HalfSpace`, `Hyperplane`, `Line2D`, and subtypes of
      `AbstractHPolygon` are preserved.
    - If `P` is an `AbstractPolytope`, then the output is an `Interval` if
      `m = 1`, an `HPolygon` if `m = 2`, and an `HPolytope` in other cases.
    - Otherwise the output is an `HPolyhedron`.

### Notes

Since the different linear-map algorithms work at the level of constraints, this
method uses dispatch on two stages: once the algorithm has been defined, first
the helper methods `_linear_map_hrep_helper` (resp. `_linear_map_vrep`) are
invoked, which dispatch on the set type. Then, each helper method calls the
concrete implementation of `_linear_map_hrep`, which dispatches on the
algorithm, and returns a list of constraints.

To simplify working with different algorithms and options, the types
`<: AbstractLinearMapAlgorithm` are used. These types are singleton type
or types that carry only the key data for the given algorithm, such as the
matrix inverse or the polyhedra backend.

New subtypes of the `AbstractPolyhedron` interface may define their own helper
methods `_linear_map_vrep` (respectively `_linear_map_hrep_helper`) for special
handling of the constraints returned by the implementations of
`_linear_map_hrep`; otherwise the fallback implementation for
`AbstractPolyhedron` is used, which instantiates an `HPolyhedron`.

### Algorithm

This method mainly implements several approaches for the linear map: inverse,
right inverse, transformation to vertex representation, variable elimination,
and variable lifting. Depending on the properties of `M` and `P`, one algorithm
may be preferable over the other. Details on the algorithms are given in the
following subsections.

Otherwise, if the algorithm argument is not specified, a default option is
chosen based on heuristics on the types and values of `M` and `P`:

- If the `"inverse"` algorithm applies, it is used.
- Otherwise, if the `"inverse_right"` algorithm applies, it is used.
- Otherwise, if the `"lift"` algorithm applies, it is used.
- Otherwise, the `"elimination"` algorithm is used.

Note that the algorithms `"inverse"` and `"inverse_right"` do not require the
external library `Polyhedra`. However, the fallback method `"elimination"`
requires `Polyhedra` as well as the library `CDDLib`.

The optional keyword arguments `inverse` and `check_invertibility` modify the
default behavior:

- If an inverse matrix is passed in `inverse`, the given algorithm is applied,
  and if none is given, either `"inverse"` or `"inverse_right"` is applied
  (in that order of preference).
- If `check_invertibility` is set to `false`, the given algorithm is applied,
  and if none is given, either `"inverse"` or `"inverse_right"` is applied
  (in that order of preference).

#### Inverse

This algorithm is invoked with the keyword argument `algorithm="inverse"`
(or `algorithm="inv"`). The algorithm requires that `M` is invertible, square,
and dense. If you know a priori that `M` is invertible, set the flag
`check_invertibility=false`, such that no extra checks are done.
Otherwise, we check the sufficient condition that the condition number
of `M` is not too high. The threshold for the condition number can be modified
from its default value, `DEFAULT_COND_TOL`, by passing a custom `cond_tol`.

The algorithm is described next. Assuming that the matrix ``M`` is invertible
(which we check via a sufficient condition,), ``y = M x`` implies
``x = \\text{inv}(M) y`` and we can transform the polyhedron
``A x ≤ b`` to the polyhedron ``A \\text{inv}(M) y ≤ b``.

If the dense condition on `M` is not satisfied, there are two suggested
workarounds: either transform to a dense matrix, i.e., calling `linear_map` with
`Matrix(M)`, or use the `"inverse_right"` algorithm, which does not compute the
inverse matrix explicitly, but uses a polyalgorithm; see the documentation
of `?\` for details.

#### Inverse-right

This algorithm is invoked with the keyword argument `algorithm="inverse_right"`
(or `algorithm="inv_right"`). This algorithm applies to square and invertible
matrices `M`. The idea is essentially the same as for the `"inverse"` algorithm;
the difference is that in `"inverse"` the full matrix inverse is computed, and
in `"inverse_right"` only the left division on the normal vectors is used. In
particular, `"inverse_right"` is good as a workaround when `M` is sparse (since
the `inv` function is not available for sparse matrices).

### Elimination

This algorithm is invoked with the keyword argument `algorithm = "elimination"`
(or `algorithm = "elim"`). The algorithm applies to any matrix `M` (invertible
or not), and any polyhedron `P` (bounded or not).

The idea is described next. If `P : Ax <= b` and `y = Mx` denote the polyhedron
and the linear map, respectively, we consider the vector `z = [y, x]`, write the
given equalities and the inequalities, and then eliminate the last x variables
(there are `length(x)` in total) using a call to `Polyhedra.eliminate` to a
backend library that can do variable elimination (typically `CDDLib` with the
`BlockElimination()` algorithm). In this way we have eliminated the "old"
variables `x` and kept the "new" or transformed variables "y".

The default elimination method is block elimination. For possible options we
refer to the documentation of Polyhedra,
[projection/elimination](https://juliapolyhedra.github.io/Polyhedra.jl/latest/projection/).

### Lift

This algorithm is invoked with the keyword argument `algorithm="lift"`.
The algorithm applies if `M` is rectangular of size `m × n` with `m > n` and
full rank (i.e., of rank `n`).

The idea is to embed the polyhedron into the `m`-dimensional space by appending
zeros, i.e. extending all constraints of `P` to `m` dimensions, and constraining
the last `m - n` dimensions to `0`. The resulting matrix is extended to an
invertible `m × m` matrix, and the algorithm using the inverse of the linear map
is applied. For technical details of extending `M` to a higher-dimensional
invertible matrix, see `ReachabilityBase.Arrays.extend`.

### Vertex representation

This algorithm is invoked with the keyword argument `algorithm="vrep"` (or
`algorithm="vrep_chull"`). If the polyhedron is bounded, the idea is to convert
it to its vertex representation and apply the linear map to each vertex.

The returned set is a polytope in vertex representation. Note that conversion of
the result back to half-space representation is not computed by default, since
this may be costly. If you use this algorithm and still want to convert back to
half-space representation, apply `tohrep` to the result of this method.
"""
function _linear_map_polyhedron(M::AbstractMatrix,
                                P::LazySet;
                                algorithm::Union{String,Nothing}=nothing,
                                check_invertibility::Bool=true,
                                cond_tol::Number=DEFAULT_COND_TOL,
                                inverse::Union{AbstractMatrix,Nothing}=nothing,
                                backend=nothing,
                                elimination_method=nothing)
    N = promote_type(eltype(M), eltype(P))
    N != eltype(P) && error("conversion between numeric types of polyhedra not " *
                            "implemented yet (see #1181)")
    if eltype(M) != eltype(P)
        M = N.(M)
    end

    got_algorithm = !isnothing(algorithm)
    got_inv = got_algorithm && (algorithm == "inv" || algorithm == "inverse")
    got_inv_right = got_algorithm && (algorithm ==
                                      "inv_right" || algorithm == "inverse_right")

    if !isnothing(inverse)
        if !got_algorithm || got_inv
            algo = LinearMapInverse(inverse)
        elseif got_inv_right
            algo = LinearMapInverseRight()
        else
            throw(ArgumentError("received an inverse matrix, but only the " *
                                "algorithms \"inverse\" and \"inverse_right\" apply, got " *
                                "$algorithm"))
        end
        return _linear_map_hrep_helper(M, P, algo)
    elseif !check_invertibility
        if got_inv || !issparse(M)
            inverse = inv(M)
            algo = LinearMapInverse(inverse)
        else
            algo = LinearMapInverseRight()
        end
        return _linear_map_hrep_helper(M, P, algo)
    end

    if !got_algorithm
        algo = _default_linear_map_algorithm(M, P; cond_tol=cond_tol)
        return _linear_map_hrep_helper(M, P, algo)

    elseif algorithm == "vrep"
        algo = LinearMapVRep(backend)
        return _linear_map_vrep(M, P, algo; apply_convex_hull=false)

    elseif algorithm == "vrep_chull"
        algo = LinearMapVRep(backend)
        return _linear_map_vrep(M, P, algo; apply_convex_hull=true)

    elseif got_inv
        check_invertibility && _check_algorithm_applies(M, P, LinearMapInverse;
                                                        cond_tol=cond_tol, throw_error=true)
        return _linear_map_hrep_helper(M, P, LinearMapInverse(inv(M)))

    elseif got_inv_right
        check_invertibility && _check_algorithm_applies(M, P, LinearMapInverseRight;
                                                        cond_tol=cond_tol, throw_error=true)
        return _linear_map_hrep_helper(M, P, LinearMapInverseRight())

    elseif algorithm == "elimination" || algorithm == "elim"
        algo = _get_elimination_instance(N, backend, elimination_method)
        return _linear_map_hrep_helper(M, P, algo)

    elseif algorithm == "lift"
        _check_algorithm_applies(M, P, LinearMapLift; throw_error=true)
        return _linear_map_hrep_helper(M, P, LinearMapLift())

    else
        throw(ArgumentError("got unknown algorithm \"$algorithm\"; available" *
                            "choices: \"inverse\", \"inverse_right\", \"lift\", " *
                            "\"elimination\", \"vrep\""))
    end
end

# TODO: merge the preconditions into _check_algorithm_applies ?
# review this method after #998
function _linear_map_vrep(M::AbstractMatrix, P::AbstractPolyhedron,
                          algo::LinearMapVRep=LinearMapVRep(nothing);
                          apply_convex_hull::Bool=false)
    require(@__MODULE__, :Polyhedra; fun_name="linear_map",
            explanation="of a $(typeof(P)) by a non-invertible matrix")

    Q = convert(VPolytope, P)
    return _linear_map_vrep(M, Q, algo; apply_convex_hull=apply_convex_hull)
end

# preconditions should have been checked in the caller function
function _linear_map_hrep(M::AbstractMatrix, P::LazySet, algo::LinearMapInverse)
    return _affine_map_inverse_hrep(algo.inverse, P)
end

# P = {y : Cy <= d}
# C(Ax + b) <= d  <=>  CAx <= d - Cb
function _affine_map_inverse_hrep(A::AbstractMatrix, P::LazySet,
                                  b::Union{AbstractVector,Nothing}=nothing)
    C_leq_d = constraints_list(P)
    constraints_res = _preallocate_constraints(C_leq_d)
    i = 1
    @inbounds for c_leq_di in C_leq_d
        cinv = vec(At_mul_B(c_leq_di.a, A))
        rhs = isnothing(b) ? c_leq_di.b : c_leq_di.b - first(At_mul_B(c_leq_di.a, b))
        if iszero(cinv)
            # constraint is redundant or infeasible
            N = eltype(cinv)
            if rhs < zero(N)
                # constraint is infeasible
                # return constraints representing empty set
                a1 = zeros(N, length(cinv))
                a1[1] = one(N)
                a2 = zeros(N, length(cinv))
                a2[1] = -one(N)
                return [HalfSpace(a1, zero(N)), HalfSpace(a2, -one(N))]
            end
        else
            constraints_res[i] = HalfSpace(cinv, rhs)
            i += 1
        end
    end
    if i <= length(constraints_res)  # there were redundant constraints, so shorten vector
        resize!(constraints_res, i - 1)
    end
    return constraints_res
end

# preconditions should have been checked in the caller function
function _linear_map_hrep(M::AbstractMatrix, P::AbstractPolyhedron,
                          algo::LinearMapInverseRight)
    constraints_P = constraints_list(P)
    constraints_MP = _preallocate_constraints(constraints_P)
    @inbounds for (i, c) in enumerate(constraints_P)
        # take left division for each constraint c, transpose(M) \ c.a
        cinv = At_ldiv_B(M, c.a)
        constraints_MP[i] = HalfSpace(cinv, c.b)
    end
    return constraints_MP
end

# preconditions should have been checked in the caller function
function _linear_map_hrep(M::AbstractMatrix, P::AbstractPolyhedron, algo::LinearMapLift)
    m, n = size(M)
    N = promote_type(eltype(M), eltype(P))

    # we extend M to an invertible m x m matrix by appending m-n columns
    # orthogonal to the column space of M
    Mext, inv_Mext = extend(M; check_rank=false)

    # append zeros to the existing constraints, in the last m-n coordinates
    # TODO: cast to common vector type instead of Vector(c.a), see #1942, #1952
    clist = constraints_list(P)
    if isempty(clist)
        cext = HalfSpace{N,Vector{N}}[]
    else
        cext = [HalfSpace(vcat(Vector(c.a), zeros(N, m - n)), c.b) for c in clist]
    end

    # now fix the last m-n coordinates to zero
    id_out = Matrix(one(N) * I, m - n, m - n)
    cext = vcat(cext, [HalfSpace(vcat(zeros(N, n), id_out[i, :]), zero(N)) for i in 1:(m - n)],
                [HalfSpace(vcat(zeros(N, n), -id_out[i, :]), zero(N)) for i in 1:(m - n)])

    Pext = HPolyhedron(cext)

    # now Mext is invertible and we can apply the inverse algorithm
    return _linear_map_hrep(Mext, Pext, LinearMapInverse(inv_Mext))
end

# If P : Ax <= b and y = Mx, we consider the vector z = [y, x], write the
# equalities and the inequalities, and then eliminate the last x variables
# (there are length(x) in total) using Polyhedra.eliminate calls
# to a backend library that can do variable elimination, typically CDDLib,
# with the BlockElimination() algorithm.
function _linear_map_hrep(M::AbstractMatrix, P::AbstractPolyhedron, algo::LinearMapElimination)
    m, n = size(M)
    N = promote_type(eltype(M), eltype(P))
    Id_neg = Matrix(-one(N) * I, m, m)
    backend = algo.backend
    method = algo.method

    # extend the polytope storing the y variables first
    # append zeros to the existing constraints, in the last m-n coordinates
    # TODO: cast to common vector type instead of hard-coding Vector(c.a), see #1942 and #1952
    clist = constraints_list(P)
    if isempty(clist)
        Ax_leq_b = Polyhedra.HalfSpace{N,Vector{N}}[]
    else
        Ax_leq_b = [Polyhedra.HalfSpace(vcat(zeros(N, m), Vector(c.a)), c.b) for c in clist]
    end
    y_eq_Mx = [Polyhedra.HyperPlane(vcat(Id_neg[i, :], Vector(M[i, :])), zero(N)) for i in 1:m]

    Phrep = Polyhedra.hrep(y_eq_Mx, Ax_leq_b)
    Phrep = polyhedron(Phrep, backend) # define concrete subtype
    Peli_block = Polyhedra.eliminate(Phrep, (m + 1):(m + n), method)
    Peli_block = Polyhedra.removeduplicates(Polyhedra.hrep(Peli_block),
                                            default_lp_solver_polyhedra(N))

    # TODO: take constraints directly -- see #1988
    return constraints_list(convert(HPolyhedron, Peli_block))
end

@inline function _preallocate_constraints(constraints::Vector{<:HalfSpace{N}}) where {N}
    return Vector{HalfSpace{N,Vector{N}}}(undef, length(constraints))
end

"""
    plot_recipe(P::AbstractPolyhedron{N}, [ε]=zero(N)) where {N}

Convert a (bounded) polyhedron to a pair `(x, y)` of points for plotting.

### Input

- `P` -- bounded polyhedron
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of points that can be plotted, where `x` is the vector of
x-coordinates and `y` is the vector of y-coordinates.

### Algorithm

We first assert that `P` is bounded (i.e., that `P` is a polytope).

One-dimensional polytopes are converted to an `Interval`.
Three-dimensional or higher-dimensional polytopes are not supported.

For two-dimensional polytopes (i.e., polygons) we compute their set of vertices
using `vertices_list` and then plot the convex hull of these vertices.
"""
function plot_recipe(P::AbstractPolyhedron{N}, ε=zero(N)) where {N}
    @assert dim(P) <= 3 "cannot plot a $(dim(P))-dimensional $(typeof(P))"
    @assert isbounded(P) "cannot plot an unbounded $(typeof(P))"

    if dim(P) == 1
        Q = convert(Interval, P)
        if diameter(Q) < _ztol(N)  # flat interval
            Q = Singleton(center(Q))
        end
        return plot_recipe(Q, ε)
    elseif dim(P) == 2
        vlist = convex_hull(vertices_list(P))
        return _plot_recipe_2d_vlist(vlist, N)
    else
        return _plot_recipe_3d_polytope(P, N)
    end
end

function _plot_recipe_2d_vlist(vlist, N)
    m = length(vlist)
    if m == 0
        @warn "received a polyhedron with no vertices during plotting"
        return plot_recipe(EmptySet{N}(2), zero(N))
    end

    x = Vector{N}(undef, m)
    y = Vector{N}(undef, m)
    @inbounds for (i, vi) in enumerate(vlist)
        x[i] = vi[1]
        y[i] = vi[2]
    end

    if m > 2
        # add first vertex to "close" the polygon
        push!(x, x[1])
        push!(y, y[1])
    end
    return x, y
end

"""
    an_element(P::AbstractPolyhedron; [solver]=default_lp_solver(eltype(P)))

Return some element of a polyhedron.

### Input

- `P`       -- polyhedron
- `solver`  -- (optional, default: `default_lp_solver(N)`) LP solver

### Output

An element of the polyhedron, or an error if the polyhedron is empty.

### Algorithm

An element is obtained by solving a feasibility linear program.
"""
function an_element(P::AbstractPolyhedron; solver=default_lp_solver(eltype(P)))
    A, b = tosimplehrep(P)
    N = eltype(P)

    lbounds, ubounds = -Inf, Inf
    sense = '<'
    obj = zeros(N, size(A, 2))
    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)

    if is_lp_optimal(lp.status)
        return lp.sol
    elseif is_lp_infeasible(lp.status)
        error("cannot return an element because the polyhedron is empty")
    else
        error("LP returned status $(lp.status) unexpectedly")
    end
end

"""
    isbounded(P::AbstractPolyhedron; [solver]=default_lp_solver(eltype(P)))

Check whether a polyhedron is bounded.

### Input

- `P`       -- polyhedron
- `solver`  -- (optional, default: `default_lp_solver(N)`) the backend used
               to solve the linear program

### Output

`true` iff the polyhedron is bounded

### Algorithm

We first check whether the polyhedron has more than `dim(P)` constraints and, if
so, whether the constraints are feasible, which is a necessary condition for
boundedness.

If so, we check boundedness via `_isbounded_stiemke`.
"""
function isbounded(P::AbstractPolyhedron; solver=default_lp_solver(eltype(P)))
    return isbounded(constraints_list(P); solver=solver)
end

function isbounded(constraints::AbstractVector{<:HalfSpace{N}};
                   solver=default_lp_solver(N)) where {N}
    if isempty(constraints)
        return false
    elseif length(constraints) <= dim(first(constraints))
        # need at least n+1 constraints to be bounded unless infeasible
        if isfeasible(constraints; solver=solver)
            return false
        end
    end
    return _isbounded_stiemke(constraints; solver=solver)
end

"""
    _isbounded_stiemke(constraints::AbstractVector{<:HalfSpace{N}};
                       solver=default_lp_solver(N),
                       check_nonempty::Bool=true) where {N}

Check whether a list of constraints is bounded using Stiemke's theorem of
alternatives.

### Input

- `constraints`    -- list of constraints
- `backend`        -- (optional, default: `default_lp_solver(N)`) the backend
                      used to solve the linear program
- `check_nonempty` -- (optional, default: `true`) if `true`, check the
                      precondition to this algorithm that `P` is non-empty

### Output

`true` iff the list of constraints is bounded.

### Notes

The list of constraints represents a polyhedron.

The algorithm calls `isempty` to check whether the polyhedron is empty.
This computation can be avoided using the `check_nonempty` flag.

### Algorithm

The algorithm is based on Stiemke's theorem of alternatives, see, e.g., [Mangasarian94](@citet).

Let the polyhedron ``P`` be given in constraint form ``Ax ≤ b``. We assume that
the polyhedron is non-empty.

Proposition 1. If ``\\ker(A)≠\\{0\\}``, then ``P`` is unbounded.

Proposition 2. Assume that ``ker(A)={0}`` and ``P`` is non-empty.
Then ``P`` is bounded if and only if the following linear program admits a
feasible solution: ``\\min∥y∥_1`` subject to ``A^Ty=0`` and ``y≥1``.
"""
function _isbounded_stiemke(constraints::AbstractVector{<:HalfSpace{N}};
                            solver=default_lp_solver(N),
                            check_nonempty::Bool=true) where {N}
    if check_nonempty && _isinfeasible(constraints; solver=solver)
        return true
    end

    A, b = tosimplehrep(constraints)
    m, n = size(A)

    if !isempty(nullspace(A))
        return false
    end

    At = copy(transpose(A))
    c = ones(N, m)
    lp = linprog(c, At, '=', zeros(n), one(N), Inf, solver)
    return is_lp_optimal(lp.status)
end

"""
    vertices_list(P::AbstractPolyhedron; check_boundedness::Bool=true)

Return the list of vertices of a polyhedron in constraint representation.

### Input

- `P`                 -- polyhedron in constraint representation
- `check_boundedness` -- (optional, default: `true`) if `true`, check whether
                         the polyhedron is bounded

### Output

The list of vertices of `P`, or an error if `P` is unbounded.

### Notes

This function throws an error if the polyhedron is unbounded. Otherwise, the
polyhedron is converted to an `HPolytope` and its list of vertices is computed.

### Examples

```jldoctest
julia> P = HPolyhedron([HalfSpace([1.0, 0.0], 1.0),
                        HalfSpace([0.0, 1.0], 1.0),
                        HalfSpace([-1.0, 0.0], 1.0),
                        HalfSpace([0.0, -1.0], 1.0)]);

julia> length(vertices_list(P))
4
```
"""
function vertices_list(P::AbstractPolyhedron; check_boundedness::Bool=true)
    if check_boundedness && !isbounded(P)
        throw(ArgumentError("the list of vertices of an unbounded " *
                            "polyhedron is not defined"))
    end
    return vertices_list(HPolytope(constraints_list(P); check_boundedness=false))
end

"""
    project(P::AbstractPolyhedron, block::AbstractVector{Int}; [kwargs...])

Concrete projection of a polyhedral set.

### Input

- `P`     -- set
- `block` -- block structure, a vector with the dimensions of interest

### Output

A polyhedron representing the projection of `P` on the dimensions specified by
`block`.
If `P` was bounded, the result is an `HPolytope`; otherwise the result is an
`HPolyhedron`.
Note that there are more specific methods for specific input types, which give a
different output type; e.g., projecting a `Ball1` results in a `Ball1`.

### Algorithm

- We first try to exploit the special case where each of the constraints of `P`
  and `block` are *compatible*, which is one of the two cases described below.
  Let `c` be a constraint of `P` and let ``D_c`` and ``D_b`` be the set of
  dimensions in which `c` resp. `block` are constrained.
  - If ``D_c ⊆ D_b``, then one can project the normal vector of `c`.
  - If ``D_c ∩ D_b = ∅``, then the constraint becomes redundant.
- In the general case, we compute the concrete linear map of the projection
  matrix associated to the given block structure.

### Examples

Consider the four-dimensional cross-polytope (unit ball in the 1-norm):

```jldoctest project_polyhedron
julia> P = convert(HPolytope, Ball1(zeros(4), 1.0));
```

All dimensions are constrained, and computing the (trivial) projection on the
whole space behaves as expected:

```jldoctest project_polyhedron
julia> constrained_dimensions(P)
1:4

julia> project(P, [1, 2, 3, 4]) == P
true
```
Each constraint of the cross polytope is constrained in all dimensions.

Now let us take a ball in the infinity norm and remove some constraints:

```jldoctest project_polyhedron
julia> B = BallInf(zeros(4), 1.0);

julia> c = constraints_list(B)[1:2]
2-element Vector{HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}}:
 HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}([1.0, 0.0, 0.0, 0.0], 1.0)
 HalfSpace{Float64, ReachabilityBase.Arrays.SingleEntryVector{Float64}}([0.0, 1.0, 0.0, 0.0], 1.0)

julia> P = HPolyhedron(c);

julia> constrained_dimensions(P)
2-element Vector{Int64}:
 1
 2
```

Finally, we take the concrete projection onto variables `1` and `2`:

```jldoctest project_polyhedron
julia> project(P, [1, 2]) |> constraints_list
2-element Vector{HalfSpace{Float64, Vector{Float64}}}:
 HalfSpace{Float64, Vector{Float64}}([1.0, 0.0], 1.0)
 HalfSpace{Float64, Vector{Float64}}([0.0, 1.0], 1.0)
```
"""
@validate function project(P::AbstractPolyhedron, block::AbstractVector{Int}; kwargs...)
    # cheap case
    clist = nothing  # allocate later
    @inbounds for c in constraints(P)
        status = _check_constrained_dimensions(c, block)
        if status == 0
            # general case
            clist = _project_polyhedron(P, block; kwargs...)
            break
        elseif status == 1
            # simple projection of half-space
            hs = HalfSpace(c.a[block], c.b)
            if isnothing(clist)
                clist = [hs]  # get the right type of the constraints
            else
                push!(clist, hs)
            end
        elseif status == -1
            # drop the constraint because it became redundant
        end
    end

    if isnothing(clist)  # set is unconstrained in the given dimensions
        N = eltype(P)
        return Universe{N}(length(block))
    end

    if isboundedtype(typeof(P))
        return HPolytope(clist; check_boundedness=false)
    else
        return HPolyhedron(clist)
    end
end

function _check_constrained_dimensions(c::HalfSpace, block)
    # 1 = constrained dimensions subset of block
    # -1 = constrained dimensions disjoint from block
    # 0 = mixed case
    status = 0
    for i in constrained_dimensions(c)
        if i in block
            # case 1: is a subset
            if status == -1
                status = 0
                break
            end
            status = 1
        else
            # case 1: is disjoint
            if status == 1
                status = 0
                break
            end
            status = -1
        end
    end
    return status
end

function _project_polyhedron(P::LazySet, block; kwargs...)
    N = eltype(P)
    M = projection_matrix(block, dim(P), N)
    πP = linear_map(M, P; kwargs...)
    return constraints_list(πP)
end

# convenience function to invert the result of `isfeasible` while still including the witness result
function _isinfeasible(constraints::AbstractVector{<:HalfSpace},
                       witness::Bool=false; solver=nothing)
    if witness
        res, w = isfeasible(constraints, witness; solver=solver)
        return !res, w
    else
        return !isfeasible(constraints, witness; solver=solver)
    end
end

"""
    addconstraint!(P::AbstractPolyhedron, constraint::HalfSpace)

Add a linear constraint to a set in constraint representation in-place.

### Input

- `P`          -- set in constraint representation
- `constraint` -- linear constraint to add

### Notes

It is left to the user to guarantee that the dimension of all linear constraints
is the same.
"""
function addconstraint!(::AbstractPolyhedron, ::HalfSpace) end

"""
    ishyperplanar(P::AbstractPolyhedron)

Determine whether a polyhedron is equivalent to a hyperplane.

### Input

- `P` -- polyhedron

### Output

`true` iff `P` is hyperplanar, i.e., consists of two linear constraints
``a·x ≤ b`` and ``-a·x ≤ -b``.
"""
function ishyperplanar(::AbstractPolyhedron) end

function extrema(P::AbstractPolyhedron)
    if dim(P) == 1
        l, h = _extrema_1d(P)
        return ([l], [h])
    end
    return _extrema_lowhigh(P)
end

@validate function extrema(P::AbstractPolyhedron, i::Int)
    if dim(P) == 1
        return _extrema_1d(P)
    end
    return _extrema_lowhigh(P, i)
end

function _extrema_1d(P::AbstractPolyhedron)
    N = eltype(P)
    l = N(-Inf)
    h = N(Inf)
    @inbounds for c in constraints(P)
        a = c.a[1]
        if a > zero(N)
            # a·x <= b  =>  x <= b/a
            h = min(h, c.b / a)
        else
            # HalfSpace `0·x <= b` is not allowed (type constraint)
            @assert a < zero(N) "invalid constraint"

            # -a·x <= -b  =>  x >= b/a
            l = max(l, c.b / a)
        end
    end
    return (l, h)
end

function low(P::AbstractPolyhedron)
    if dim(P) == 1
        l = _low_1d(P)
        return [l]
    end
    return _low(P)
end

@validate function low(P::AbstractPolyhedron, i::Int)
    if dim(P) == 1
        return _low_1d(P)
    end
    return _low(P, i)
end

function _low_1d(P::AbstractPolyhedron)
    N = eltype(P)
    l = N(-Inf)
    @inbounds for c in constraints(P)
        a = c.a[1]
        if a < zero(N)
            # -a·x <= -b  =>  x >= b/a
            l = max(l, c.b / a)
        end
    end
    return l
end

function high(P::AbstractPolyhedron)
    if dim(P) == 1
        h = _high_1d(P)
        return [h]
    end
    return _high(P)
end

@validate function high(P::AbstractPolyhedron, i::Int)
    if dim(P) == 1
        return _high_1d(P)
    end
    return _high(P, i)
end

function _high_1d(P::AbstractPolyhedron)
    N = eltype(P)
    h = N(Inf)
    @inbounds for c in constraints(P)
        a = c.a[1]
        if a > zero(N)
            # a·x <= b  =>  x <= b/a
            h = min(h, c.b / a)
        end
    end
    return h
end

# create n+1 contradicting constraints
# n times ⋀_i x_i ≤ 0
# 1 times ∑_i x_i ≥ 1
# Note: constraints are sorted CCW in 2D
function _infeasible_constraints_list(n::Int; N=Float64)
    c_sum = HalfSpace(fill(N(-1), n), N(-1))  # ∑_i x_i ≥ 1
    clist = Vector{typeof(c_sum)}(undef, n + 1)
    @inbounds for i in 1:n
        a = zeros(N, n)
        a[i] = one(N)
        clist[i] = HalfSpace(a, N(0))  # x_i ≤ 0
    end
    @inbounds clist[n + 1] = c_sum
    return clist
end
