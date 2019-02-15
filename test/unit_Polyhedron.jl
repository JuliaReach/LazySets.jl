using SparseArrays
using LazySets: isinvertible

global test_suite_polyhedra

for N in [Float64, Rational{Int}, Float32]
    # random polyhedron
    rand(HPolyhedron)

    # constructor from matrix and vector
    A = [N(1) N(2); N(-1) N(1)]
    b = [N(1), N(2)]
    p = HPolyhedron(A, b)
    c = p.constraints
    @test c isa Vector{LinearConstraint{N}}
    @test c[1].a == N[1, 2] && c[1].b == N(1)
    @test c[2].a == N[-1, 1] && c[2].b == N(2)

    # convert back to matrix and vector
    A2, b2 = tosimplehrep(p)
    @test A == A2 && b == b2

    # 2D HPolyhedron
    p = HPolyhedron{N}()
    c1 = LinearConstraint(N[2, 2], N(12))
    c2 = LinearConstraint(N[-3, 3], N(6))
    c3 = LinearConstraint(N[-1, -1], N(0))
    c4 = LinearConstraint(N[2, -4], N(0))
    addconstraint!(p, c3)
    addconstraint!(p, c1)
    addconstraint!(p, c4)
    addconstraint!(p, c2)

    # tosimplehrep for other polytopic types
    A, b = tosimplehrep(BallInf(N[0], N(2)))
    @test (A == hcat(N[1; -1]) || A == hcat(N[-1; 1])) && b == N[2, 2]
    # tosimplehrep from list of constraints
    A, b = tosimplehrep([HalfSpace(N[1], N(2)), HalfSpace(N[-1], N(2))])
    @test (A == hcat(N[1; -1]) || A == hcat(N[-1; 1])) && b == N[2, 2]

    # support vector
    d = N[1, 0]
    @test σ(d, p) == N[4, 2]
    d = N[0, 1]
    @test σ(d, p) == N[2, 4]
    d = N[-1, 0]
    @test σ(d, p) == N[-1, 1]
    d = N[0, -1]
    @test σ(d, p) == N[0, 0]

    # support vector of polyhedron with no constraints
    @test σ(N[1], HPolyhedron{N}()) == N[Inf]

    # boundedness
    @test isbounded(p)
    @test !isbounded(HPolyhedron{N}())

    # membership
    @test ∈(N[5 / 4, 7 / 4], p)
    @test !∈(N[4, 1], p)

    # constrained dimensions
    @test constrained_dimensions(p) == [1, 2]
    @test constrained_dimensions(
        HPolyhedron{N}([LinearConstraint(N[1, 0], N(1))])) == [1]

    # concrete linear map with invertible matrix
    linear_map(N[2 3; 1 2], p)

    if test_suite_polyhedra
        # conversion to and from Polyhedra's VRep data structure
        cl = constraints_list(HPolyhedron(polyhedron(p)))
        @test length(p.constraints) == length(cl)

        # convert hyperrectangle to a HPolyhedron
        H = Hyperrectangle(N[1, 1], N[2, 2])
        P = convert(HPolyhedron, H)
        @test_throws ArgumentError vertices_list(P) # the vertices list is not defined, see #820

        # checking for emptiness
        P = HPolyhedron([LinearConstraint(N[1, 0], N(0))])    # x <= 0
        @test !isempty(P)
        addconstraint!(P, LinearConstraint(N[-1, 0], N(-1)))  # x >= 1
        @test isempty(P)

        # concrete linear map with noninvertible matrix throws an error
        @test_throws ArgumentError linear_map(N[2 3; 0 0], P)
    end
end

# default Float64 constructors
unconstrained_HPolyhedron = HPolyhedron()
@test unconstrained_HPolyhedron isa HPolyhedron{Float64}

# Polyhedra tests that only work with Float64
if test_suite_polyhedra
    for N in [Float64]
        # support function/vector
        d = N[1, 0]
        p_unbounded = HPolyhedron([LinearConstraint(N[-1, 0], N(0))])
        @test σ(d, p_unbounded) == N[Inf, 0]
        @test ρ(d, p_unbounded) == N(Inf)
        p_infeasible = HPolyhedron([LinearConstraint(N[1], N(0)),
                                  LinearConstraint(N[-1], N(-1))])
        @test_throws ErrorException σ(N[1], p_infeasible)

        # intersection
        # TODO these polyhedra are empty. do the tests make any sense?
        A = [N(0) N(1); N(1) N(0); N(2) N(2)]
        b = N[0, 0, 1]
        p1 = HPolyhedron(A, b)
        A = [N(0) N(-1); N(-1) N(0); N(1) N(1)]
        b = N[-0.25, -0.25, 0]
        p2 = HPolyhedron(A, b)
        cap = intersection(p1, p2)
        # currently broken, see #565

        # intersection with polytope
        A = [N(1) N(0); N(0) N(1); N(-1) N(-1)]
        b = N[1, 1, -1]
        p = HPolyhedron(A, b)
        b = BallInf(N[1, 1], N(0.5))
        cap = intersection(p, b)
        cap = intersection(b, p)

        # intersection with half-space
        hs = HalfSpace(N[2, 2], N(-1))
        c1 = constraints_list(intersection(p1, hs))
        c2 = constraints_list(intersection(hs, p1))
        @test length(c1) == 3 && ispermutation(c1, c2)

        # convex hull
        ch = convex_hull(p1, p2)
        # currently broken, see #566

        # Cartesian product
        A = [N(1) N(-1)]'
        b = N[1, 0]
        p1 = HPolyhedron(A, b)
        p2 = HPolyhedron(A, b)
        cp = cartesian_product(p1, p2)
        cl = constraints_list(cp)
        @test length(cl) == 4

        # vertices_list
        A = [N(1) N(-1)]'
        b = N[1, 0]
        p = HPolyhedron(A, b)
        @test_throws ArgumentError vertices_list(p)

        # tovrep from HPolyhedron
        A = [N(0) N(-1); N(-1) N(0); N(1) N(1)]
        b = N[-0.25, -0.25, -0]
        P = HPolyhedron(A, b)
        @test tohrep(P) isa HPolyhedron # test no-op

        # removing a redundant constraint in a 2D polyhedron (see #565)
        A = [0. 1.;    # y <= 0
             1. 0.;    # x <= 0
             2. 2.]    # 2x + 2y <= 1   (this constraint is implied by the other two)
        b = [0., 0., 1.];
        p1 = HPolyhedron(A, b);
        remove_redundant_constraints!(p1)
        Ar, br = tosimplehrep(p1)
        @test Ar == A[1:2, :] && br == b[1:2]

        # removing redundant constraints from a list of constraints
        constraints = [HalfSpace(N[1], N(1)), HalfSpace(N[1], N(0))]
        constraints2 = remove_redundant_constraints(constraints)
        result = remove_redundant_constraints!(constraints)
        @test result
        @test ispermutation(constraints, constraints2)
        constraints = [HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-1)), HalfSpace(N[-1], N(-1))]
        result = remove_redundant_constraints!(constraints)
        @test !result

        # checking for empty intersection (also test symmetric calls)
        P = convert(HPolytope, BallInf(zeros(N, 2), N(1)))
        Q = convert(HPolytope, BallInf(ones(N, 2), N(1)))
        R = HPolyhedron([HalfSpace(N[1, 0], N(3)), HalfSpace(N[-1, 0], N(-2))])
        res, w = isdisjoint(P, Q, true)
        @test !isdisjoint(P, Q) && !res && w ∈ P && w ∈ Q
        res, w = isdisjoint(Q, P, true)
        @test !isdisjoint(Q, P) && !res && w ∈ P && w ∈ Q
        res, w = isdisjoint(Q, R, true)
        @test !isdisjoint(Q, R) && !res && w ∈ Q && w ∈ R
        res, w = isdisjoint(R, Q, true)
        @test !isdisjoint(R, Q) && !res && w ∈ Q && w ∈ R
        res, w = isdisjoint(P, R, true)
        @test isdisjoint(P, R) && res && w == N[]
        res, w = isdisjoint(R, P, true)
        @test isdisjoint(R, P) && res && w == N[]

        Punbdd = HPolyhedron([HalfSpace([0.68, 1.22], -0.76),
                              HalfSpace{Float64}([-0.75, -0.46], 0.68)])
        @assert !isbounded(Punbdd)

        Pbdd = HPolyhedron([HalfSpace([0.68, 1.22], -0.76),
                            HalfSpace{Float64}([-0.75, -0.46], 0.68),
                            HalfSpace{Float64}([1.72, 0.33], 0.37),
                            HalfSpace{Float64}([-1.60, -0.41], 0.67),
                            HalfSpace{Float64}([-0.44, 0.06], 0.78)])
        @assert isbounded(Pbdd)

        Mnotinv = [1.0 0.0; 2.0 0.0]
        @assert !isinvertible(Mnotinv)

        Minv = [1.0 2.0; -1.0 0.4]
        @assert isinvertible(Minv)

        # invertible matrix times a bounded polyhedron 
        L = linear_map(Minv, Pbdd)
        @test L isa HPolyhedron

        # invertible matrix times an unbounded polyhedron 
        L = linear_map(Minv, Punbdd)
        @test L isa HPolyhedron

        # not invertible matrix times a bounded polyhedron 
        L = linear_map(Mnotinv, Pbdd) # Requires Polyhedra because it works on vertices
        @test L isa VPolytope

        # not invertible matrix times an unbounded polyhedron
        @test_throws ArgumentError linear_map(Mnotinv, Punbdd)

        # check that we can use sparse matrices as well ; Requires SparseArrays
        L = linear_map(sparse(Minv), Pbdd)
        @test L isa HPolyhedron
        L = linear_map(sparse(Minv), Punbdd)
        @test L isa HPolyhedron
        L = linear_map(sparse(Mnotinv), Pbdd) # Requires Polyhedra because it works on vertices
        @test L isa VPolytope
        @test_throws ArgumentError linear_map(sparse(Mnotinv), Punbdd)
    end
end
