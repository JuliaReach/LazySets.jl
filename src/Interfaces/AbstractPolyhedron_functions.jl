import Base.∈

export constrained_dimensions,
       tosimplehrep,
       remove_redundant_constraints,
       remove_redundant_constraints!,
       linear_map,
       chebyshev_center,
       an_element,
       vertices_list,
       singleton_list

# default LP solver for floating-point numbers
function default_lp_solver(N::Type{<:AbstractFloat})
    GLPKSolverLP(method=:Simplex)
end

# default LP solver for rational numbers
function default_lp_solver(N::Type{<:Rational})
    GLPKSolverLP(method=:Exact)
end

# fallback method
function default_polyhedra_backend(P, N)
    require(:Polyhedra; fun_name="default_polyhedra_backend")
    error("no default backend for numeric type $N")
end

"""
    ∈(x::AbstractVector{N}, P::AbstractPolyhedron{N}) where {N<:Real}

Check whether a given point is contained in a polyhedron.

### Input

- `x` -- point/vector
- `P` -- polyhedron

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation checks if the point lies inside each defining half-space.
"""
function ∈(x::AbstractVector{N}, P::AbstractPolyhedron{N}) where {N<:Real}
    @assert length(x) == dim(P) "a $(length(x))-dimensional point cannot be " *
        "an element of a $(dim(P))-dimensional set"

    for c in constraints_list(P)
        if dot(c.a, x) > c.b
            return false
        end
    end
    return true
end

"""
    isuniversal(P::AbstractPolyhedron{N}, [witness]::Bool=false
               ) where {N<:Real}

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
function isuniversal(P::AbstractPolyhedron{N}, witness::Bool=false
                    ) where {N<:Real}
    constraints = constraints_list(P)
    if isempty(constraints)
        return witness ? (true, N[]) : true
    else
        return witness ? isuniversal(constraints[1], true) : false
    end
end

"""
    constrained_dimensions(P::AbstractPolyhedron) where {N<:Real}

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
    tosimplehrep(constraints::AbstractVector{LC})
        where {N<:Real, LC<:LinearConstraint{N}}

Return the simple H-representation ``Ax ≤ b`` from a list of linear constraints.

### Input

- `constraints` -- a list of linear constraints

### Output

The tuple `(A, b)` where `A` is the matrix of normal directions and `b` is the
vector of offsets.
"""
function tosimplehrep(constraints::AbstractVector{LC}
                     ) where {N<:Real, LC<:LinearConstraint{N}}
    n = length(constraints)
    if n == 0
        A = Matrix{N}(undef, 0, 0)
        b = Vector{N}(undef, 0)
        return (A, b)
    end
    A = zeros(N, n, dim(first(constraints)))
    b = zeros(N, n)
    @inbounds begin
        for (i, Pi) in enumerate(constraints)
            A[i, :] = Pi.a
            b[i] = Pi.b
        end
    end
    return (A, b)
end

"""
     remove_redundant_constraints!(
         constraints::AbstractVector{<:LinearConstraint{N}};
         [backend]=default_lp_solver(N)) where {N<:Real}

Remove the redundant constraints of a given list of linear constraints; the list
is updated in-place.

### Input

- `constraints` -- list of constraints
- `backend`     -- (optional, default: `default_lp_solver(N)`) the backend used
                   to solve the linear program
### Output

`true` if the function was successful and the list of constraints `constraints`
is modified by removing the redundant constraints, and `false` only if the
constraints are infeasible.

### Notes

Note that the result may be `true` even if the constraints are infeasible.
For example, ``x ≤ 0 && x ≥ 1`` will return `true` without removing any
constraint.
To check if the constraints are infeasible, use
`isempty(HPolyhedron(constraints)`.

### Algorithm

If there are `m` constraints in `n` dimensions, this function checks one by one
if each of the `m` constraints is implied by the remaining ones.

To check if the `k`-th constraint is redundant, an LP is formulated using the
constraints that have not yet being removed.
If, at an intermediate step, it is detected that a subgroup of the constraints
is infeasible, this function returns `false`.
If the calculation finished successfully, this function returns `true`.

For details, see [Fukuda's Polyhedra
FAQ](https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node24.html).
"""
function remove_redundant_constraints!(
        constraints::AbstractVector{<:LinearConstraint{N}};
        backend=default_lp_solver(N)) where {N<:Real}

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
        if lp.status == :Infeasible
            # the polyhedron is empty
            return false
        elseif lp.status == :Optimal
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

    deleteat!(constraints, setdiff(1:m, non_redundant_indices))
    return true
end

"""
    remove_redundant_constraints(
        constraints::AbstractVector{<:LinearConstraint{N}};
        backend=default_lp_solver(N)) where {N<:Real}

Remove the redundant constraints of a given list of linear constraints.

### Input

- `constraints` -- list of constraints
- `backend`     -- (optional, default: `default_lp_solver(N)`) the backend used
                   to solve the linear program
### Output

The list of constraints with the redundant ones removed, or an empty set if the
constraints are infeasible.

### Algorithm

See `remove_redundant_constraints!(::AbstractVector{<:LinearConstraint})` for
details.
"""
function remove_redundant_constraints(
        constraints::AbstractVector{<:LinearConstraint{N}};
        backend=default_lp_solver(N)) where {N<:Real}
    constraints_copy = copy(constraints)
    if remove_redundant_constraints!(constraints_copy, backend=backend)
        return constraints_copy
    else  # the constraints are infeasible
        return EmptySet{N}(dim(constraints[1]))
    end
end

# =======================
# Concrete linear map
# =======================

struct LinearMapInverse{T, MT<:AbstractMatrix{T}} <: AbstractLinearMapAlgorithm
    inverse::MT
end

struct LinearMapInverseRight <: AbstractLinearMapAlgorithm
    #
end

struct LinearMapLift <: AbstractLinearMapAlgorithm
    #
end

struct LinearMapElimination{T, S} <: AbstractLinearMapAlgorithm
    backend::T
    method::S
end

struct LinearMapVRep <: AbstractLinearMapAlgorithm
    #
end

function _check_algorithm_applies(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                ::Type{LinearMapInverse}; cond_tol=DEFAULT_COND_TOL, throw_error=false) where {N}

    inv_condition = issquare(M) && isinvertible(M; cond_tol=cond_tol)
    if !inv_condition
        throw_error && throw(ArgumentError("algorithm \"inverse\" requires an " *
                            "invertible matrix"))
        return false
    end

    dense_condition = !issparse(M)
    if !dense_condition
        throw_error && throw(ArgumentError("the inverse of a sparse matrix is not " *
            "available; either convert your matrix to a dense matrix with `Matrix(M)`, " *
            "or try the \"inverse_right\" algorithm"))
        return false
    end
    return true
end

function _check_algorithm_applies(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
            ::Type{LinearMapInverseRight}; cond_tol=DEFAULT_COND_TOL, throw_error=false) where {N}
    inv_condition = issquare(M) && isinvertible(M; cond_tol=cond_tol)
    if !inv_condition
        throw_error && throw(ArgumentError("algorithm \"inverse_right\" requires an " *
                            "invertible matrix"))
        return false
    end
    return true
end

function _check_algorithm_applies(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                                   ::Type{LinearMapLift}; throw_error=false) where {N}

    m, n = size(M)
    size_condition = m > n
    if !size_condition
        throw_error && throw(ArgumentError("algorithm \"lift\" requires that the number " *
            "of rows of the linear map is greater than the number of columns, but " *
            "they are of size $m and $n respectively"))
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

function _check_algorithm_applies(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                                   ::Type{LinearMapVRep}; throw_error=false) where {N}

    # TODO: the second check should be !isbounded(P) but this may be expensive;
    # see also #998 and #1926
    is_polytopic = P isa AbstractPolytope
    if !is_polytopic
        throw_error && throw(ArgumentError("algorithm \"vrep\" requires that the " *
            "list of vertices of the polyhedron is available, but it is not for a polyhedron " *
            "of type $(typeof(P))"))
        return false
    end
    return true
end

function _get_elimination_instance(N, backend, elimination_method)
    require(:Polyhedra; fun_name="linear_map with elimination")
    if backend == nothing
        require(:CDDLib; fun_name="linear_map with elimination")
        backend = default_cddlib_backend(N)
    end
    if elimination_method == nothing
        elimination_method = Polyhedra.BlockElimination()
    end
    return LinearMapElimination(backend, elimination_method)
end

function _default_linear_map_algorithm(M::AbstractMatrix{N}, P::AbstractPolyhedron{N};
                cond_tol=DEFAULT_COND_TOL, backend=nothing, elimination_method=nothing) where {N}

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
    linear_map(M::AbstractMatrix{N},
               P::AbstractPolyhedron{N};
               [algorithm]::Union{String, Nothing}=nothing,
               [check_invertibility]::Bool=true,
               [cond_tol]::Number=DEFAULT_COND_TOL,
               [inverse]::Union{AbstractMatrix{N}, Nothing}=nothing,
               [backend]=nothing,
               [elimination_method]=nothing) where {N<:Real}

Concrete linear map of a polyhedral set.

### Input

- `M`         -- matrix
- `P`         -- polyhedral set
- `algorithm` -- (optional; default: `nothing`) algorithm to be used; for the
                 description see the Algorithm section below; possible choices are:

    - `"inverse"`, alias: `"inv"`
    - `"inverse_right"`, alias: `"inv_right"`
    - `"elimination"`, alias: `"elim"`
    - `"lift"`
    - `"vrep"`

- `check_invertibility` -- (optional, default: `true`) if `true` check whether
                           given matrix `M` is invertible; set to `false` only
                           if you know that `M` is invertible
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
    - Otherwise, the output is an `Interval` if `m = 1`, a `VPolygon` if `m = 2`
      and a `VPolytope` in other cases.

- If the invertibility criterion was used:

    - The types of `HalfSpace`, `Hyperplane`, `Line` and `AbstractHPolygon` are
      preserved.
    - If `P` is an `AbstractPolytope`, then the output is an `Interval` if `m = 1`,
      an `HPolygon` if `m = 2` and an `HPolytope` in other cases.
    - Otherwise, the output is an `HPolyhedron`.

### Notes

Since the different linear map algorithms work at the level of constraints (not sets
representations), this function uses dispatch on two stages: once the algorithm
has been defined, first the helper functions `_linear_map_hrep_helper` (resp.
`_linear_map_vrep`) are invoked, which dispatch on the set type. Then, each helper
function calls the concrete implementation of `_linear_map_hrep`, which dispatches
on the algorithm, and returns a list of constraints.

To simplify working with different algorithms and options, the
types `<: AbstractLinearMapAlgorithm` are used. These types are singleton type
or types that carry only the key data for the given algorithm, such as the matrix
inverse or the polyhedra backend.

New subtypes of the `AbstractPolyhedron` interface may define their own helper functions
`_linear_map_vrep`, respectively `_linear_map_hrep_helper` for special handling
of the constraints returned by the implementations of `_linear_map_hrep`;
otherwise the fallback implementation for `AbstractPolyhedron` is used, which
instantiates an `HPolyhedron`.

### Algorithm

This function mainly implements several approaches for the linear map: inverse,
right inverse, transformation to the vertex representation, variable elimination,
and variable lifting. Depending on the properties of `M` and `P`, one algorithm
may be preferable over the other. Details on the algorithms are given in the
following subsections.

Otherwise, if the algorithm argument is not specified, a default option is chosen based
on heuristics on the types and values of `M` and `P`:

- If the `"inverse"` algorithm applies, it is used.
- If the `"inverse_right"` algorithm applies, it is used.
- Otherwise, if the `"lift"` algorithm applies, it is used.
- Otherwise, the `"elimination"` algorithm is used.

Note that `"inverse"` does not require the external library `Polyhedra`, and neither
does `"inverse_right"`. However, the fallback method `"elimination"` requires
`Polyhedra` as well as the library `CDDLib`.

The optional keyword arguments `inverse` and `check_invertibility`
modify the default behavior:

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
`check_invertibility=false`, such that no extra checks are done within `linear_map`.
Otherwise, we check the sufficient condition that the condition number
of `M` is not too high. The threshold for the condition number can be modified
from its default value, `DEFAULT_COND_TOL`, by passing a custom `cond_tol`.

The algorithm is described next. Assuming that the matrix ``M`` is invertible
(which we check via a sufficient condition,), ``y = M x`` implies
``x = \\text{inv}(M) y`` and we can transform the polyhedron
``A x ≤ b`` to the polyhedron ``A \\text{inv}(M) y ≤ b``.

If the dense condition on `M` is not fullfilled, there are two suggested
workarounds: either transform to dense matrix, i.e. calling `linear_map` with
`Matrix(M)`, or use the `"inverse_right"` algorithm, which does not compute the
inverse matrix explicitly, but uses a polyalgorithm; see the documentation
of `?\` for details.

#### Inverse-right

This algorithm is invoked with the keyword argument `algorithm="inverse_right"`
(or `algorithm="inv_right"`). This algorithm applies to square and invertible matrices
`M`. The idea is essentially the same as for the `"inverse"` algorithm; the difference
is that in `"inverse"` the full matrix inverse is computed, and in `"inverse_right"`
only the left division on the normal vectors is used. In particular, `"inverse_right"`
is good as a workaround when `M` is sparse (since the `inv` function is not available
for sparse matrices).

### Elimination

This algorithm is invoked with the keyword argument `algorithm = "elimination"` or
`algorithm = "elim"`. The algorithm applies to any matrix `M` (invertible or not),
and any polyhedron `P` (bounded or not).

The idea is described next. If `P : Ax <= b` and `y = Mx` denote the polyhedron
and the linear map respectively, we consider the vector `z = [y, x]`, write the
given equalities and the inequalities, and then eliminate the last x variables
(there are `length(x)` in total) using a call to `Polyhedra.eliminate` to a backend
library that can do variable elimination, typically `CDDLib` with the
`BlockElimination()` algorithm. In this way we have eliminated the "old" variables
`x` and kept the "new" or transformed variables "y".

The default elimination method is block elimination. For possible options we refer
to the documentation of Polyhedra,
[projection/elimination](https://juliapolyhedra.github.io/Polyhedra.jl/latest/projection/).

### Lift

This algorithm is invoked with the keyword argument `algorithm="lift"`.
The algorithm applies if `M` is rectangular of size `m × n` with `m > n` and
full rank (i.e. of rank `n`).

The idea is to embed the polyhedron into the `m`-dimensional space by appending zeros,
i.e. extending all constraints of `P` to `m` dimensions, and constraining the last
`m - n` dimensions to `0`. The matrix resulting matrix is extended to an invertible
`m × m` matrix and the algorithm using the inverse of the linear map is applied.
For the technical details of the extension of `M` to a higher-dimensional
invertible matrix, see `LazySets.Arrays.extend`.

### Vertex representation

This algorithm is invoked with the keyword argument `algorithm="vrep"`.
The idea is to convert the polyhedron to its vertex representation and apply the
linear map to each vertex of `P`.

The returned set is a polytope in vertex representation. Note that conversion of
the result back to half-space representation is not computed by default, since this
may be costly. If you used this algorithm and still want to convert back to
half-space representation, apply `tohrep` to the result of this function.
Note that this method only works for bounded polyhedra.
"""
function linear_map(M::AbstractMatrix{N},
                    P::AbstractPolyhedron{N};
                    algorithm::Union{String, Nothing}=nothing,
                    check_invertibility::Bool=true,
                    cond_tol::Number=DEFAULT_COND_TOL,
                    inverse::Union{AbstractMatrix{N}, Nothing}=nothing,
                    backend=nothing,
                    elimination_method=nothing) where {N<:Real}

   size(M, 2) != dim(P) && throw(ArgumentError("a linear map of size $(size(M)) " *
                            "cannot be applied to a set of dimension $(dim(P))"))

    got_algorithm = algorithm != nothing
    got_inv = got_algorithm && (algorithm == "inv" || algorithm == "inverse")
    got_inv_right = got_algorithm && (algorithm == "inv_right" || algorithm == "inverse_right")

    if inverse != nothing
        if !got_algorithm || got_inv
            algo = LinearMapInverse(inverse)
        elseif got_inv_right
            algo = LinearMapInverseRight()
        else
            throw(ArgumentError("received an inverse matrix but only the algorithms " *
                                "\"inverse\" and \"inverse_right\" apply, got $algorithm"))
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
        return _linear_map_vrep(M, P)

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
        throw(ArgumentError("got unknown algorithm \"$algorithm\"; " *
            "available choices: \"inverse\", \"inverse_right\", \"lift\", " *
            "\"elimination\", \"vrep\""))
    end
end

# handle different numeric types
function linear_map(M::AbstractMatrix{NM},
                    P::AbstractPolyhedron{NP};
                    kwargs...) where {NM<:Real, NP<:Real}
    N = promote_type(NM, NP)
    if N != NP
        error("conversion between numeric types of polyhedra not implemented " *
            "yet (see #1181)")
    else
        return linear_map(N.(M), P; kwargs...)
    end
end

# TODO: merge the preconditions into _check_algorithm_applies ?
# review this method after #998
function _linear_map_vrep(M::AbstractMatrix{N}, P::AbstractPolyhedron{N}) where {N<:Real}
    if !isbounded(P)
        throw(ArgumentError("the linear map in vertex representation for an " *
            "unbounded set is not defined"))
    end
    require(:Polyhedra; fun_name="linear_map",
            explanation="of a $(typeof(P)) by a non-invertible matrix")
    # since P is bounded, we pass an HPolytope and then convert it to vertex representation
    P = tovrep(HPolytope(constraints_list(P), check_boundedness=false))
    return _linear_map_vrep(M, P)
end

# generic function for the AbstractPolyhedron interface => returns an HPolyhedron
function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                                 algo::AbstractLinearMapAlgorithm) where {N<:Real}
    constraints = _linear_map_hrep(M, P, algo)
    return HPolyhedron(constraints)
end

# preconditions should have been checked in the caller function
function _linear_map_hrep(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                          algo::LinearMapInverse) where {N}
    inverse = algo.inverse
    constraints_P = constraints_list(P)
    constraints_MP = _preallocate_constraints(constraints_P)
    @inbounds for (i, c) in enumerate(constraints_P)
        cinv = vec(_At_mul_B(c.a, inverse))
        constraints_MP[i] = LinearConstraint(cinv, c.b)
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
        cinv = _At_ldiv_B(M, c.a)
        constraints_MP[i] = LinearConstraint(cinv, c.b)
    end
    return constraints_MP
end

# preconditions should have been checked in the caller function
function _linear_map_hrep(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                          algo::LinearMapLift) where {N<:Real}
    m, n = size(M)

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
function _linear_map_hrep(M::AbstractMatrix{N}, P::AbstractPolyhedron{N},
                          algo::LinearMapElimination) where {N<:Real}
    m, n = size(M)
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
    Peli_block = Polyhedra.removeduplicates(Polyhedra.hrep(Peli_block))

    # TODO: take constraints directly -- see #1988
    return constraints_list(HPolyhedron(Peli_block))
end

@inline function _preallocate_constraints(constraints::Vector{LinearConstraint{N, VN}}) where {N, VN<:AbstractVector{N}}
    return Vector{LinearConstraint{N, Vector{N}}}(undef, length(constraints))
end

"""
    plot_recipe(P::AbstractPolyhedron{N}, [ε]::N=zero(N)) where {N<:Real}

Convert a (bounded) polyhedron to a pair `(x, y)` of points for plotting.

### Input

- `P` -- bounded polyhedron
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of points that can be plotted.

### Algorithm

We first assert that `P` is bounded (i.e., that `P` is a polytope).

One-dimensional polytopes are converted to an `Interval`.
Three-dimensional or higher-dimensional polytopes are not supported.

For two-dimensional polytopes (i.e., polygons) we compute their set of vertices
using `vertices_list` and then plot the convex hull of these vertices.
"""
function plot_recipe(P::AbstractPolyhedron{N}, ε::N=zero(N)) where {N<:Real}
    @assert dim(P) <= 2 "cannot plot a $(dim(P))-dimensional $(typeof(P))"
    @assert isbounded(P) "cannot plot an unbounded $(typeof(P))"

    if dim(P) == 1
        return plot_recipe(convert(Interval, P), ε)
    else
        vlist = transpose(hcat(convex_hull(vertices_list(P))...))
        if isempty(vlist)
            @warn "received a polyhedron with no vertices during plotting"
            return plot_recipe(EmptySet{N}(2), ε)
        end
        x, y = vlist[:, 1], vlist[:, 2]

        if length(x) > 1
            # add first vertex to "close" the polygon
            push!(x, x[1])
            push!(y, y[1])
        end

        return x, y
    end
end

"""
    chebyshev_center(P::AbstractPolyhedron{N};
                     [get_radius]::Bool=false,
                     [backend]=default_polyhedra_backend(P, N),
                     [solver]=JuMP.with_optimizer(GLPK.Optimizer)
                     ) where {N<:AbstractFloat}

Compute the [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
of a polytope.

### Input

- `P`          -- polytope
- `get_radius` -- (optional; default: `false`) option to additionally return the
                  radius of the largest ball enclosed by `P` around the
                  Chebyshev center
- `backend`    -- (optional; default: `default_polyhedra_backend(P, N)`) the
                  backend for polyhedral computations
- `solver`     -- (optional; default: `JuMP.with_optimizer(GLPK.Optimizer)`) the
                  LP solver passed to `Polyhedra`

### Output

If `get_radius` is `false`, the result is the Chebyshev center of `P`.
If `get_radius` is `true`, the result is the pair `(c, r)` where `c` is the
Chebyshev center of `P` and `r` is the radius of the largest ball with center
`c` enclosed by `P`.

### Notes

The Chebyshev center is the center of a largest Euclidean ball enclosed by `P`.
In general, the center of such a ball is not unique (but the radius is).
"""
function chebyshev_center(P::AbstractPolyhedron{N};
                          get_radius::Bool=false,
                          backend=default_polyhedra_backend(P, N),
                          solver=JuMP.with_optimizer(GLPK.Optimizer)
                         ) where {N<:AbstractFloat}
    require(:Polyhedra; fun_name="chebyshev_center")
    # convert to HPolyhedron to ensure `polyhedron` is applicable (see #1505)
    Q = polyhedron(convert(HPolyhedron, P); backend=backend)
    c, r = Polyhedra.chebyshevcenter(Q, solver)

    if get_radius
        return c, r
    end
    return c
end

"""
    an_element(P::AbstractPolyhedron{N}) where {N<:Real}

Return some element of a convex set.

### Input

- `P` -- polyhedron

### Output

An element of a polyhedron.

### Algorithm

An element of the polyhedron is obtained by evaluating its support vector along
direction ``[1, 0, …, 0]``.
"""
function an_element(P::AbstractPolyhedron{N}) where {N<:Real}
    n = dim(P)
    if n == -1
        throw(ArgumentError("the dimension of this polyhedron is not defined, " *
                            "hence `an_element` is not available"))
    end
    e₁ = SingleEntryVector(1, n, one(N))
    return σ(e₁, P)
end

"""
    vertices_list(P::AbstractPolyhedron; check_boundedness::Bool=true)

Return the list of vertices of a polyhedron in constraint representation.

### Input

- `P`                 -- polyhedron in constraint representation
- `check_boundedness` -- (optional, default: `true`) if `true`, check whether the
                         polyhedron is bounded

### Output

The list of vertices of `P`, or an error if `P` is unbounded.

### Notes

This function returns an error if the polyhedron is unbounded. Otherwise,
the polyhedron is converted to an `HPolytope` and its list of vertices is computed.

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
    return vertices_list(HPolytope(constraints_list(P), check_boundedness=false))
end

"""
    singleton_list(P::AbstractPolyhedron; check_boundedness::Bool=true)

Return the vertices of a polyhedron in H-representation as a list of singletons.

### Input

- `P`                 -- polyhedron in constraint representation
- `check_boundedness` -- (optional, default: `true`) if `true`, check whether the
                         polyhedron is bounded

### Output

The list of vertices of `P`, as `Singleton`, or an error if `P` is unbounded.

### Notes

This function returns an error if the polyhedron is unbounded. Otherwise,
the polyhedron is converted to an `HPolytpe` and its list of vertices is computed.
"""
function singleton_list(P::AbstractPolyhedron; check_boundedness::Bool=true)
    return [Singleton(x) for x in vertices_list(P)]
end
