import Base.∈

export constrained_dimensions,
       tosimplehrep,
       remove_redundant_constraints,
       remove_redundant_constraints!,
       linear_map,
       an_element,
       vertices_list

# default LP solver for floating-point numbers
function default_lp_solver(N::Type{<:AbstractFloat})
    JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(method=GLPK.SIMPLEX))
end

# default LP solver for rational numbers
function default_lp_solver(N::Type{<:Rational})
    JuMP.optimizer_with_attributes(() -> GLPK.Optimizer(method=GLPK.EXACT))
end

# helper function given two possibly different numeric types
function default_lp_solver(M::Type{<:Number}, N::Type{<:Number})
    return default_lp_solver(promote_type(M, N))
end

# Polyhedra backend (fallback method)
function default_polyhedra_backend(P::LazySet{N}) where {N}
    require(@__MODULE__, :Polyhedra; fun_name="default_polyhedra_backend")
    error("no default backend for numeric type $N")
end

# default LP solver for Polyhedra (fallback method)
# NOTE: exists in parallel to `default_lp_solver` because we use different
# interfaces (see #1493)
function default_lp_solver_polyhedra(N, varargs...)
    require(@__MODULE__, :Polyhedra; fun_name="default_lp_solver_polyhedra")
    error("no default solver for numeric type $N")
end

isconvextype(::Type{<:AbstractPolyhedron}) = true

is_polyhedral(::AbstractPolyhedron) = true

"""
    ∈(x::AbstractVector, P::AbstractPolyhedron)

Check whether a given point is contained in a polyhedron.

### Input

- `x` -- point/vector
- `P` -- polyhedron

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation checks if the point lies inside each defining half-space.
"""
function ∈(x::AbstractVector, P::AbstractPolyhedron)
    @assert length(x) == dim(P) "a $(length(x))-dimensional point cannot be " *
        "an element of a $(dim(P))-dimensional set"

    for c in constraints_list(P)
        if !_leq(dot(c.a, x), c.b)
            return false
        end
    end
    return true
end

"""
    isuniversal(P::AbstractPolyhedron{N}, [witness]::Bool=false) where {N}

Check whether a polyhedron is universal.

### Input

- `P`       -- polyhedron
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``P`` is universal
* If `witness` option is activated:
  * `(true, [])` iff ``P`` is universal
  * `(false, v)` iff ``P`` is not universal and ``v ∉ P``

### Algorithm

`P` is universal iff it has no constraints.

A witness is produced using `isuniversal(H)` where `H` is the first linear
constraint of `P`.
"""
function isuniversal(P::AbstractPolyhedron{N}, witness::Bool=false) where {N}
    c = constraints(P)
    if isempty(c)
        return witness ? (true, N[]) : true
    else
        return witness ? isuniversal(first(c), true) : false
    end
end

"""
    constrained_dimensions(P::AbstractPolyhedron)

Return the indices in which a polyhedron is constrained.

### Input

- `P` -- polyhedron

### Output

A vector of ascending indices `i` such that the polyhedron is constrained in
dimension `i`.

### Examples

A 2D polyhedron with constraint ``x1 ≥ 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(P::AbstractPolyhedron)
    constraints = constraints_list(P)
    if isempty(constraints)
        return Int[]
    end
    zero_indices = zeros(Int, dim(P))
    for constraint in constraints
        for i in constrained_dimensions(constraint)
            zero_indices[i] = i
        end
    end
    return filter(x -> x != 0, zero_indices)
end

"""
    tosimplehrep(constraints::AbstractVector{LC};
                 [n]::Int=0) where {N, LC<:HalfSpace{N}}

Return the simple H-representation ``Ax ≤ b`` from a list of linear constraints.

### Input

- `constraints` -- a list of linear constraints
- `n`           -- (optional; default: `0`) dimension of the constraints

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` is the
vector of offsets.

### Notes

The parameter `n` can be used to create a matrix with no constraints but a
non-zero dimension.
"""
function tosimplehrep(constraints::AbstractVector{LC};
                      n::Int=0) where {N, LC<:HalfSpace{N}}
    m = length(constraints)
    if m == 0
        A = Matrix{N}(undef, 0, n)
        b = Vector{N}(undef, 0)
        return (A, b)
    end
    if n <= 0
        n = dim(first(constraints))
    end
    A = zeros(N, m, n)
    b = zeros(N, m)
    @inbounds begin
        for (i, Pi) in enumerate(constraints)
            A[i, :] = Pi.a
            b[i] = Pi.b
        end
    end
    return (A, b)
end

"""
    remove_redundant_constraints!(constraints::AbstractVector{S};
                                  [backend]=nothing) where {S<:HalfSpace}

Remove the redundant constraints of a given list of linear constraints; the list
is updated in-place.

### Input

- `constraints` -- list of constraints
- `backend`     -- (optional, default: `nothing`) the backend used to solve the
                   linear program
### Output

`true` if the removal was successful and the list of constraints `constraints`
is modified by removing the redundant constraints, and `false` only if the
constraints are infeasible.

### Notes

Note that the result may be `true` even if the constraints are infeasible.
For example, ``x ≤ 0 && x ≥ 1`` will return `true` without removing any
constraint.
To check if the constraints are infeasible, use
`isempty(HPolyhedron(constraints))`.

If `backend` is `nothing`, it defaults to `default_lp_solver(N)`.

### Algorithm

If there are `m` constraints in `n` dimensions, this function checks one by one
if each of the `m` constraints is implied by the remaining ones.

To check if the `k`-th constraint is redundant, an LP is formulated using the
constraints that have not yet been removed.
If, at an intermediate step, it is detected that a subgroup of the constraints
is infeasible, this function returns `false`.
If the calculation finished successfully, this function returns `true`.

For details, see [Fukuda's Polyhedra
FAQ](https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node24.html).
"""
function remove_redundant_constraints!(constraints::AbstractVector{S};
                                       backend=nothing) where {S<:HalfSpace}
    if isempty(constraints)
        return true
    end
    N = eltype(first(constraints))
    if isnothing(backend)
        backend = default_lp_solver(N)
    end
    A, b = tosimplehrep(constraints)
    m, n = size(A)
    non_redundant_indices = 1:m

    i = 1 # counter over reduced constraints

    for j in 1:m    # loop over original constraints
        α = A[j, :]
        Ar = A[non_redundant_indices, :]
        br = b[non_redundant_indices]
        br[i] = b[j] + one(N)
        lp = linprog(-α, Ar, '<', br, -Inf, Inf, backend)
        if is_lp_infeasible(lp.status)
            # the polyhedron is empty
            return false
        elseif is_lp_optimal(lp.status)
            objval = -lp.objval
            if _leq(objval, b[j])
                # the constraint is redundant
                non_redundant_indices = setdiff(non_redundant_indices, j)
            else
                # the constraint is not redundant
                i = i+1
            end
        else
            error("LP is not optimal; the status of the LP is $(lp.status)")
        end
    end

    idx_delete = setdiff(1:m, non_redundant_indices)
    if !isempty(idx_delete)
        deleteat!(constraints, idx_delete)
        return true
    else
        return false
    end
end

"""
    remove_redundant_constraints(constraints::AbstractVector{S};
                                 backend=nothing) where {S<:HalfSpace}

Remove the redundant constraints of a given list of linear constraints.

### Input

- `constraints` -- list of constraints
- `backend`     -- (optional, default: `nothing`) the backend used to solve the
                   linear program
### Output

The list of constraints with the redundant ones removed, or an empty set if the
constraints are infeasible.

### Notes

If `backend` is `nothing`, it defaults to `default_lp_solver(N)`.

### Algorithm

See `remove_redundant_constraints!(::AbstractVector{<:HalfSpace})` for
details.
"""
function remove_redundant_constraints(constraints::AbstractVector{S};
                                      backend=nothing) where {S<:HalfSpace}
    constraints_copy = copy(constraints)
    if remove_redundant_constraints!(constraints_copy, backend=backend)
        return constraints_copy
    else  # the constraints are infeasible
        N = eltype(first(constraints))
        return EmptySet{N}(dim(constraints[1]))
    end
end

struct LinearMapInverse{T, MT<:AbstractMatrix{T}} <: AbstractLinearMapAlgorithm
    inverse::MT
end

struct LinearMapInverseRight <: AbstractLinearMapAlgorithm end

struct LinearMapLift <: AbstractLinearMapAlgorithm end

struct LinearMapElimination{T, S} <: AbstractLinearMapAlgorithm
    backend::T
    method::S
end

struct LinearMapVRep{T} <: AbstractLinearMapAlgorithm
    backend::T
end

function _check_algorithm_applies(M::AbstractMatrix{N},
                                  P::AbstractPolyhedron{N},
                                  ::Type{LinearMapInverse};
                                  cond_tol=DEFAULT_COND_TOL,
                                  throw_error=false) where {N}

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

function _check_algorithm_applies(M::AbstractMatrix{N},
                                  P::AbstractPolyhedron{N},
                                  ::Type{LinearMapInverseRight};
                                  cond_tol=DEFAULT_COND_TOL,
                                  throw_error=false) where {N}
    inv_condition = issquare(M) && isinvertible(M; cond_tol=cond_tol)
    if !inv_condition
        throw_error && throw(ArgumentError("algorithm \"inverse_right\" " *
                                           "requires an invertible matrix"))
        return false
    end
    return true
end

function _check_algorithm_applies(M::AbstractMatrix{N},
                                  P::AbstractPolyhedron{N},
                                  ::Type{LinearMapLift};
                                  throw_error=false) where {N}

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

function _check_algorithm_applies(M::AbstractMatrix{N},
                                  P::AbstractPolyhedron{N},
                                  ::Type{LinearMapVRep};
                                  throw_error=false) where {N}

    if !isboundedtype(typeof(P))
        throw_error && throw(ArgumentError("algorithm \"vrep\" requires a " *
            "polytope and its list of vertices, but the set is a $(typeof(P))"))
        return false
    end
    return true
end

function _get_elimination_instance(N, backend, elimination_method)
    require(@__MODULE__, :Polyhedra; fun_name="linear_map with elimination")
    if backend == nothing
        require(@__MODULE__, :CDDLib; fun_name="linear_map with elimination")
        backend = default_cddlib_backend(N)
    end
    if elimination_method == nothing
        elimination_method = Polyhedra.BlockElimination()
    end
    return LinearMapElimination(backend, elimination_method)
end

function _default_linear_map_algorithm(M::AbstractMatrix{N},
                                       P::AbstractPolyhedron{N};
                                       cond_tol=DEFAULT_COND_TOL,
                                       backend=nothing,
                                       elimination_method=nothing) where {N}

    if _check_algorithm_applies(M, P, LinearMapInverse, cond_tol=cond_tol)
        algo = LinearMapInverse(inv(M))
    elseif _check_algorithm_applies(M, P, LinearMapLift)
        algo = LinearMapLift()
    else
        algo = _get_elimination_instance(N, backend, elimination_method)
    end
    return algo
end

"""
    linear_map(M::AbstractMatrix{NM},
               P::AbstractPolyhedron{NP};
               [algorithm]::Union{String, Nothing}=nothing,
               [check_invertibility]::Bool=true,
               [cond_tol]::Number=DEFAULT_COND_TOL,
               [inverse]::Union{AbstractMatrix{N}, Nothing}=nothing,
               [backend]=nothing,
               [elimination_method]=nothing) where {NM, NP}

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

If the dense condition on `M` is not fullfilled, there are two suggested
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
function linear_map(M::AbstractMatrix{NM},
                    P::AbstractPolyhedron{NP};
                    algorithm::Union{String, Nothing}=nothing,
                    check_invertibility::Bool=true,
                    cond_tol::Number=DEFAULT_COND_TOL,
                    inverse::Union{AbstractMatrix, Nothing}=nothing,
                    backend=nothing,
                    elimination_method=nothing) where {NM, NP}
    N = promote_type(NM, NP)
    N != NP && error("conversion between numeric types of polyhedra not " *
        "implemented yet (see #1181)")
    if NM != NP
        M = N.(M)
    end

    size(M, 2) != dim(P) && throw(ArgumentError("a linear map of size " *
                "$(size(M)) cannot be applied to a set of dimension $(dim(P))"))

    got_algorithm = algorithm != nothing
    got_inv = got_algorithm && (algorithm == "inv" || algorithm == "inverse")
    got_inv_right = got_algorithm && (algorithm ==
        "inv_right" || algorithm == "inverse_right")

    if inverse != nothing
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
        _check_algorithm_applies(M, P, LinearMapLift, throw_error=true)
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
    if !isbounded(P)
        throw(ArgumentError("the linear map in vertex representation for an " *
                            "unbounded set is not defined"))
    end
    require(@__MODULE__, :Polyhedra; fun_name="linear_map",
            explanation="of a $(typeof(P)) by a non-invertible matrix")
    # since P is bounded, we pass an HPolytope and then convert it to vertex
    # representation

    P_hpoly = HPolytope(constraints_list(P), check_boundedness=false)
    backend = algo.backend
    if backend == nothing
        backend = default_polyhedra_backend(P)
    end
    P = tovrep(P_hpoly, backend=backend)
    return _linear_map_vrep(M, P, algo; apply_convex_hull=apply_convex_hull)
end

# generic function for the AbstractPolyhedron interface => returns an HPolyhedron
function _linear_map_hrep_helper(M::AbstractMatrix, P::AbstractPolyhedron,
                                 algo::AbstractLinearMapAlgorithm)
    constraints = _linear_map_hrep(M, P, algo)
    return HPolyhedron(constraints)
end

# preconditions should have been checked in the caller function
function _linear_map_hrep(M::AbstractMatrix, P::AbstractPolyhedron,
                          algo::LinearMapInverse)
    return _linear_map_inverse_hrep(algo.inverse, P)
end

function linear_map_inverse(Minv::AbstractMatrix, P::AbstractPolyhedron)
    @assert size(Minv, 1) == dim(P) "a linear map of size $(size(Minv)) " *
        "cannot be applied to a set of dimension $(dim(P))"
    constraints = _linear_map_inverse_hrep(Minv, P)
    return HPolyhedron(constraints)
end

function _linear_map_inverse_hrep(Minv::AbstractMatrix, P::AbstractPolyhedron)
    constraints_P = constraints_list(P)
    constraints_MP = _preallocate_constraints(constraints_P)
    @inbounds for (i, c) in enumerate(constraints_P)
        cinv = vec(At_mul_B(c.a, Minv))
        constraints_MP[i] = HalfSpace(cinv, c.b)
    end
    return constraints_MP
end

# preconditions should have been checked in the caller function
function _linear_map_hrep(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                          algo::LinearMapInverseRight) where {N}
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
function _linear_map_hrep(M::AbstractMatrix{NM}, P::AbstractPolyhedron{NP},
                          algo::LinearMapLift) where {NM, NP}
    m, n = size(M)
    N = promote_type(NM, NP)

    # we extend M to an invertible m x m matrix by appending m-n columns
    # orthogonal to the column space of M
    Mext, inv_Mext = extend(M, check_rank=false)

    # append zeros to the existing constraints, in the last m-n coordinates
    # TODO: cast to common vector type instead of Vector(c.a), see #1942, #1952
    cext = [HalfSpace(vcat(Vector(c.a), zeros(N, m-n)), c.b) for c in constraints_list(P)]

    # now fix the last m-n coordinates to zero
    id_out = Matrix(one(N)*I, m-n, m-n)
    cext = vcat(cext, [HalfSpace(vcat(zeros(N, n), id_out[i, :]), zero(N)) for i in 1:(m-n)],
                      [HalfSpace(vcat(zeros(N, n), -id_out[i, :]), zero(N)) for i in 1:(m-n)])

    Pext = HPolyhedron(cext)

    # now Mext is invertible and we can apply the inverse algorithm
    return _linear_map_hrep(Mext, Pext, LinearMapInverse(inv_Mext))
end

# If P : Ax <= b and y = Mx, we consider the vector z = [y, x], write the
# equalities and the inequalities, and then eliminate the last x variables
# (there are length(x) in total) using Polyhedra.eliminate calls
# to a backend library that can do variable elimination, typically CDDLib,
# with the BlockElimination() algorithm.
function _linear_map_hrep(M::AbstractMatrix{NM}, P::AbstractPolyhedron{NP},
                          algo::LinearMapElimination) where {NM, NP}
    m, n = size(M)
    N = promote_type(NM, NP)
    ₋Id_m = Matrix(-one(N)*I, m, m)
    backend = algo.backend
    method = algo.method

    # extend the polytope storing the y variables first
    # append zeros to the existing constraints, in the last m-n coordinates
    # TODO: cast to common vector type instead of hard-coding Vector(c.a), see #1942 and #1952
    Ax_leq_b = [Polyhedra.HalfSpace(vcat(zeros(N, m), Vector(c.a)), c.b) for c in constraints_list(P)]
    y_eq_Mx = [Polyhedra.HyperPlane(vcat(₋Id_m[i, :], Vector(M[i, :])), zero(N)) for i in 1:m]

    Phrep = Polyhedra.hrep(y_eq_Mx, Ax_leq_b)
    Phrep = polyhedron(Phrep, backend) # define concrete subtype
    Peli_block = Polyhedra.eliminate(Phrep, (m+1):(m+n), method)
    Peli_block = Polyhedra.removeduplicates(Polyhedra.hrep(Peli_block),
                                            default_lp_solver_polyhedra(N))

    # TODO: take constraints directly -- see #1988
    return constraints_list(HPolyhedron(Peli_block))
end

@inline function _preallocate_constraints(constraints::Vector{<:HalfSpace{N}}) where {N}
    return Vector{HalfSpace{N, Vector{N}}}(undef, length(constraints))
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
    @assert dim(P) <= 2 "cannot plot a $(dim(P))-dimensional $(typeof(P))"
    @assert isbounded(P) "cannot plot an unbounded $(typeof(P))"

    if dim(P) == 1
        Q = convert(Interval, P)
        if diameter(Q) < _ztol(N)  # flat interval
            Q = Singleton(center(Q))
        end
        return plot_recipe(Q, ε)
    else
        vlist = convex_hull(vertices_list(P))
        return _plot_recipe_2d_vlist(vlist, N)
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
    an_element(P::AbstractPolyhedron{N};
               [solver]=default_lp_solver(N)) where {N}

Return some element of a polyhedron.

### Input

- `P`       -- polyhedron
- `solver`  -- (optional, default: `default_lp_solver(N)`) LP solver

### Output

An element of the polyhedron, or an error if the polyhedron is empty.

### Algorithm

An element is obtained by solving a feasibility linear program.
"""
function an_element(P::AbstractPolyhedron{N};
                    solver=default_lp_solver(N)) where {N}

    A, b = tosimplehrep(P)

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
    isbounded(P::AbstractPolyhedron{N}; [solver]=default_lp_solver(N)) where {N}

Check whether a polyhedron is bounded.

### Input

- `P`       -- polyhedron
- `solver`  -- (optional, default: `default_lp_solver(N)`) the backend used
               to solve the linear program

### Output

`true` iff the polyhedron is bounded

### Algorithm

We first check if the polyhedron has more than `dim(P)` constraints, which is a
necessary condition for boundedness.

If so, we check boundedness via `_isbounded_stiemke`.
"""
function isbounded(P::AbstractPolyhedron{N}; solver=default_lp_solver(N)) where {N}
    return isbounded(constraints_list(P); solver=solver)
end

function isbounded(constraints::AbstractVector{<:HalfSpace{N}};
                   solver=default_lp_solver(N)) where {N}
    if isempty(constraints)
        return false
    elseif length(constraints) <= dim(first(constraints))
        return false  # need at least n+1 constraints to be bounded
    end
    return _isbounded_stiemke(constraints, solver=solver)
end

"""
    _isbounded_stiemke(constraints::AbstractVector{<:HalfSpace{N}};
                       solver=LazySets.default_lp_solver(N),
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

The algorithm is based on Stiemke's theorem of alternatives, see, e.g., [1].

Let the polyhedron ``P`` be given in constraint form ``Ax ≤ b``. We assume that
the polyhedron is non-empty.

Proposition 1. If ``\\ker(A)≠\\{0\\}``, then ``P`` is unbounded.

Proposition 2. Assume that ``ker(A)={0}`` and ``P`` is non-empty.
Then ``P`` is bounded if and only if the following linear program admits a
feasible solution: ``\\min∥y∥_1`` subject to ``A^Ty=0`` and ``y≥1``.

[1] Mangasarian, Olvi L. *Nonlinear programming.*
    Society for Industrial and Applied Mathematics, 1994.
"""
function _isbounded_stiemke(constraints::AbstractVector{<:HalfSpace{N}};
                            solver=LazySets.default_lp_solver(N),
                            check_nonempty::Bool=true) where {N}
    if check_nonempty && isempty(HPolyhedron(constraints))
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
    if check_boundedness && !isboundedtype(typeof(P)) && !isbounded(P)
        throw(ArgumentError("the list of vertices of an unbounded " *
                            "polyhedron is not defined"))
    end
    return vertices_list(HPolytope(constraints_list(P), check_boundedness=false))
end

"""
    project(P::AbstractPolyhedron{N}, block::AbstractVector{Int};
            [kwargs...]) where {N}

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
4-element Vector{Int64}:
 1
 2
 3
 4

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
function project(P::AbstractPolyhedron{N}, block::AbstractVector{Int};
                 kwargs...) where {N}
    general_case = false

    # cheap case
    clist = nothing  # allocate later
    @inbounds for c in constraints(P)
        status = _check_constrained_dimensions(c, block)
        if status == 0
            general_case = true
            break
        elseif status == 1
            # simple projection of half-space
            hs = HalfSpace(c.a[block], c.b)
            if clist == nothing
                clist = [hs]  # get the right type of the constraints
            else
                push!(clist, hs)
            end
        elseif status == -1
            # drop the constraint because it became redundant
        end
    end

    # general case
    if general_case
        clist = _project_polyhedron(P, block; kwargs...)
    end

    if isnothing(clist)  # set is unconstrained in the given dimensions
        return Universe{N}(dim(P))
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

function _project_polyhedron(P::LazySet{N}, block; kwargs...) where {N}
    M = projection_matrix(block, dim(P), N)
    πP = linear_map(M, P; kwargs...)
    return constraints_list(πP)
end
