using LazySets: AbstractPolyhedron, LinearMapElimination, constraints_list,
                default_cddlib_backend, polyhedron,
                remove_redundant_constraints!
using LazySets.HPolyhedronModule: HPolyhedron
using Polyhedra: EliminationAlgorithm, FourierMotzkin, HRep, convexhull,
                 eliminate, hrep, removeduplicates, removehredundancy!
using ReachabilityBase.Basetype: basetype
using ReachabilityBase.Subtypes: subtypes
import Base: convert
import LazySets: default_lp_solver_polyhedra, _linear_map_hrep

function convert(::Type{HPolyhedron}, P::HRep)
    return _convert_HPoly(HPolyhedron, P)
end

# This internal helper function computes the Minkowski sum of two polyhedra in
# simple H-representation, P = {x : Ax <= b} and Q = {x : Cx <= d}, using
# projection methods. See the documentation of `minkowski_sum` for details.
function _minkowski_sum_hrep(A::AbstractMatrix, b::AbstractVector,
                             C::AbstractMatrix, d::AbstractVector;
                             backend=nothing, algorithm=nothing, prune=true)
    if isnothing(backend)
        N = promote_type(eltype(A), eltype(b), eltype(C), eltype(d))
        backend = default_cddlib_backend(N)
    end

    if isnothing(algorithm)
        algorithm = FourierMotzkin()
    elseif !(algorithm isa EliminationAlgorithm)  # NOTE: this is an internal function
        throw(ArgumentError("algorithm $algorithm is not a valid elimination algorithm; " *
                            "choose among any of $(subtypes(EliminationAlgorithm))"))
    end

    mP, nP = size(A)
    mQ, nQ = size(C)
    E = [zeros(N, mP, nQ) A; C -C]
    f = [b; d]
    PQ = HPolyhedron(E, f)
    PQ_cdd = polyhedron(PQ; backend=backend)
    W_cdd = eliminate(PQ_cdd, (nP + 1):(2nP), algorithm)
    W = convert(HPolyhedron, W_cdd)
    if prune
        success = remove_redundant_constraints!(W)
        if !success
            throw(ArgumentError("the constraints corresponding to the Minkowski sum of the " *
                                "given sets are infeasible"))
        end
    end
    return W
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

    Phrep = hrep(y_eq_Mx, Ax_leq_b)
    Phrep = Polyhedra.polyhedron(Phrep, backend) # define concrete subtype
    Peli_block = eliminate(Phrep, (m + 1):(m + n), method)
    Peli_block = removeduplicates(hrep(Peli_block), default_lp_solver_polyhedra(N))

    # TODO: take constraints directly -- see #1988
    return constraints_list(convert(HPolyhedron, Peli_block))
end
