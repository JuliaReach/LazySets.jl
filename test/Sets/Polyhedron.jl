global test_suite_polyhedra

using LazySets: _isbounded_stiemke, _isbounded_unit_dimensions

for N in [Float64, Rational{Int}, Float32]
    # random polyhedron
    rand(HPolyhedron)

    p_univ = HPolyhedron{N}()

    # constructor from matrix and vector
    A = [N(1) N(2); N(-1) N(1)]
    b = [N(1), N(2)]
    p = HPolyhedron(A, b)
    c = p.constraints
    @test c isa Vector{LinearConstraint{N, Vector{N}}}
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

    # support vector of polyhedron with no constraints
    @test σ(N[1], p_univ) == N[Inf]
    @test σ(N[1], p_univ; algorithm="halfspace_direction") == N[Inf]

    # is_polyhedral
    @test is_polyhedral(p)

    # universality
    @test !isuniversal(p)
    res, w = isuniversal(p, true)
    @test !res && w ∉ p
    @test isuniversal(p_univ) && isuniversal(p_univ, true) == (true, N[])

    # membership
    @test N[5 / 4, 7 / 4] ∈ p && N[4, 1] ∉ p

    # constrained dimensions
    @test constrained_dimensions(p) == [1, 2]
    @test constrained_dimensions(
        HPolyhedron([LinearConstraint(N[1, 0], N(1))])) == [1]
    @test constrained_dimensions(HPolyhedron()) == Int[]

    # concrete linear map with invertible matrix
    linear_map(N[2 3; 1 2], p)

    # translation
    p2 = translate(p, N[1, 2])
    @test p2 isa HPolyhedron{N} && ispermutation(constraints_list(p2),
        [HalfSpace(N[2, 2], N(18)), HalfSpace(N[-3, 3], N(9)),
         HalfSpace(N[-1, -1], N(-3)), HalfSpace(N[2, -4], N(-6))])

    # constraints iterator
    @test ispermutation(collect(constraints(p)), constraints_list(p))

    # vertices iterator
    @test ispermutation(collect(vertices(p)), vertices_list(p))

    # equivalence check
    p1 = HPolyhedron([HalfSpace(N[1], N(1))])
    p2 = HPolyhedron([HalfSpace(N[1], N(-1))])
    @test isequivalent(p1, p1) && !isequivalent(p1, p2)

    if test_suite_polyhedra
        # conversion to and from Polyhedra's VRep data structure
        cl = constraints_list(HPolyhedron(polyhedron(p)))
        @test length(p.constraints) == length(cl)

        # convert hyperrectangle to a HPolyhedron
        H = Hyperrectangle(N[1, 1], N[2, 2])
        P = convert(HPolyhedron, H)
        @test length(vertices_list(P)) == 4

        # checking for emptiness
        P = HPolyhedron([LinearConstraint(N[1, 0], N(0))])    # x <= 0
        @test !isempty(P)

        # concrete linear map with non-invertible matrix
        if N == Float64
            @test linear_map(N[2 3; 0 0], P) isa HPolyhedron
        end
    end

    if !test_suite_polyhedra
        # concrete linear map of a bounded polyhedron by a non-invertible matrix
        # throws an assertion error, since tovrep(HPolytope(...)) is required
        H = Hyperrectangle(N[1, 1], N[2, 2])
        P = convert(HPolyhedron, H)
        if N != Float32
            # conversion to vrep with Float32 fails
            @test_throws ArgumentError linear_map(N[2 3; 0 0], P, algorithm="vrep")
        end

        if N != Rational{Int} # in floating-point we can use elimination
            lm = linear_map(N[2 3; 0 0], P, algorithm="elimination")
            @test lm isa HPolyhedron{Float64}

            B = N[4e8 2; 0 1]
            P = CartesianProduct(BallInf(N[0.01], N(0.08)), Singleton(N[1.0]))
            lm = linear_map(B, P)
            @test lm isa HPolytope{Float64}
        end
    end

    if test_suite_polyhedra
        # robustness of empty set (see issue #2532)
        s1 = HalfSpace(N[-1.0], -N(1.0000000000000002)) # x  >= 1.0000000000000002
        s2 = HalfSpace(N[1.0], N(1.0)) # x <= 1
        P = s1 ∩ s2

        if N == Float64
            # can't prove emptiness on Float64
            @test isempty(P) == false

        elseif N == Rational{Int}
            # can prove emptiness using exact arithmetic
            @test isempty(P, use_polyhedra_interface=true, backend=CDDLib.Library(:exact)) == true
        end
    end
end

# default Float64 constructors
@test HPolyhedron() isa HPolyhedron{Float64}

# tests that only work with Float64 and Float32
for N in [Float64, Float32]
    # normalization
    p1 = HPolyhedron([HalfSpace(N[1e5], N(3e5)), HalfSpace(N[-2e5], N(4e5))])
    p2 = normalize(p1)
    for hs in constraints_list(p2)
        @test norm(hs.a) == N(1)
    end

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

    p_univ = HPolyhedron{N}()

    # boundedness
    @test !isbounded(p_univ) && !isboundedtype(typeof(p_univ))
    @test isbounded(p)
    @test !isbounded(HPolyhedron([HalfSpace(N[1, 0], N(1))]))

    @test _isbounded_stiemke(constraints_list(p_univ))
    @test _isbounded_stiemke(constraints_list(p))
    @test !_isbounded_stiemke([HalfSpace(N[1, 0], N(1))])

    @test _isbounded_unit_dimensions(p_univ)
    @test _isbounded_unit_dimensions(p)
    @test !_isbounded_unit_dimensions(HPolyhedron([HalfSpace(N[1, 0], N(1))]))
end

# Polyhedra tests that only work with Float64
for N in [Float64]
    p = HPolyhedron{N}()
    c1 = LinearConstraint(N[2, 2], N(12))
    c2 = LinearConstraint(N[-3, 3], N(6))
    c3 = LinearConstraint(N[-1, -1], N(0))
    c4 = LinearConstraint(N[2, -4], N(0))
    addconstraint!(p, c3)
    addconstraint!(p, c1)
    addconstraint!(p, c4)
    addconstraint!(p, c2)

    # support vector
    d = N[1, 0]
    @test σ(d, p) == N[4, 2]
    @test σ(d, p; algorithm="halfspace_direction") == N[4, 2]
    @test ρ(d, p) == 4.0
    @test ρ(d, p; algorithm="halfspace_direction") == 4.0
    d = N[0, 1]
    @test σ(d, p) == N[2, 4]
    @test σ(d, p; algorithm="halfspace_direction") == N[2, 4]
    @test ρ(d, p) == 4.0
    @test ρ(d, p; algorithm="halfspace_direction") == 4.0
    d = N[-1, 0]
    @test σ(d, p) == N[-1, 1]
    @test σ(d, p; algorithm="halfspace_direction") == N[-1, 1]
    @test ρ(d, p) == 1.0
    @test ρ(d, p; algorithm="halfspace_direction") == 1.0
    d = N[0, -1]
    @test σ(d, p) == N[0, 0]
    @test σ(d, p; algorithm="halfspace_direction") == N[0, 0]
    @test ρ(d, p) == 0.0
    @test ρ(d, p; algorithm="halfspace_direction") == 0.0

    # membership
    @test [Inf, Inf] ∉ p

    # an_element
    P = HPolyhedron([HalfSpace(N[3//50, -4//10], N(1)),
                     HalfSpace(N[-1//50, 1//10], N(-1))])
    @test an_element(P) ∈ P

    # an_element for an unbounded polyhedron
    P = HPolyhedron([HalfSpace(N[-1, 0], N(-1))])
    y = an_element(P)
    # check that all entries are finite
    @test all(!isinf, y)
    # check that the points belong to P
    @test y ∈ P

    # boundedness
    @test isbounded(p)

    if test_suite_polyhedra
        p_unbounded = HPolyhedron([LinearConstraint(N[-1, 0], N(0))])
        p_infeasible = HPolyhedron([LinearConstraint(N[1], N(0)),
                                    LinearConstraint(N[-1], N(-1))])

        # support function/vector
        d = N[1, 0]
        @test σ(d, p_unbounded)[1] == N(Inf)
        @test ρ(d, p_unbounded) == N(Inf)
        @test σ(d, p_unbounded; algorithm="halfspace_direction")[1] == N(Inf)
        @test ρ(d, p_unbounded; algorithm="halfspace_direction") == N(Inf)

        d = N[-1, 1]
        @test σ(d, p_unbounded)[2] == N(Inf)
        @test ρ(d, p_unbounded) == N(Inf)
        @test σ(d, p_unbounded; algorithm="halfspace_direction")[2] == N(Inf)
        @test ρ(d, p_unbounded; algorithm="halfspace_direction") == N(Inf)

        @test_throws ErrorException σ(N[1], p_infeasible)

        # isempty
        @test !isempty(p_unbounded)
        @test isempty(p_infeasible)
        P = HPolyhedron([LinearConstraint(N[1, 0], N(0))])    # x <= 0
        addconstraint!(P, LinearConstraint(N[-1, 0], N(-1)))  # x >= 1
        @test isempty(P)

        # intersection
        A = N[0 1; 1 0]
        b = N[1, 1]
        p1 = HPolyhedron(A, b)
        A = N[0 -1; -1 0]
        b = N[0, 0]
        p2 = HPolyhedron(A, b)
        cap = intersection(p1, p2)
        @test ispermutation(constraints_list(cap),
                            constraints_list(BallInf(N[0.5, 0.5], N(0.5))))

        # intersection with polytope
        A = N[1 0;    # x <= 1
              0 1;    # y <= 1
              -1 -1]  # x + y >= 1
        b = N[1, 1, -1]
        p = HPolyhedron(A, b)
        b = BallInf(N[1, 1], N(0.5))
        cap1 = intersection(p, b)
        cap2 = intersection(b, p)
        @test ispermutation(constraints_list(cap1), constraints_list(cap2)) &&
              ispermutation(constraints_list(cap1),
                            constraints_list(BallInf(N[0.75, 0.75], N(0.25))))

        # intersection with half-space
        hs = HalfSpace(N[2, 2], N(-1))
        c1 = constraints_list(intersection(p1, hs))
        c2 = constraints_list(intersection(hs, p1))
        @test length(c1) == 3 && ispermutation(c1, c2)

        # convex hull
        ch = convex_hull(p1, p2)
        @test ch isa HPolyhedron{N} && isempty(constraints_list(ch))

        # Cartesian product
        A = N[1 -1]'
        b = N[1, 0]
        p1 = HPolyhedron(A, b)
        p2 = copy(p1)
        cp = cartesian_product(p1, p2)
        cl = constraints_list(cp)
        @test length(cl) == 4

        # vertices_list
        vertices_list(p1) ≈ [N[1.0], N[0.0]]

        # tovrep from HPolyhedron
        @test tohrep(p1) isa HPolyhedron{N} # test no-op

        # removing a redundant constraint in a 2D polyhedron (see #565)
        A = N[0 1;  # y <= 0
              1 0;  # x <= 0
              2 2]  # 2x + 2y <= 1 (the constraint is implied by the other two)
        b = N[0, 0, 1]
        p1 = HPolyhedron(A, b)
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
        # removing redundant constraints does not alter the original polyhedron
        P = HPolyhedron([HalfSpace(N[1], N(0)), HalfSpace(N[1], N(1))])
        Q = remove_redundant_constraints(P)
        @test length(P.constraints) == 2 && length(Q.constraints) == 1

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

        Punbdd = HPolyhedron([HalfSpace(N[0.68, 1.22], N(-0.76)),
                              HalfSpace(N[-0.75, -0.46], N(0.68))])
        @assert !isbounded(Punbdd)

        Pbdd = HPolyhedron([HalfSpace(N[0.68, 1.22], N(-0.76)),
                            HalfSpace(N[-0.75, -0.46], N(0.68)),
                            HalfSpace(N[1.72, 0.33], N(0.37)),
                            HalfSpace(N[-1.60, -0.41], N(0.67)),
                            HalfSpace(N[-0.44, 0.06], N(0.78))])
        @assert isbounded(Pbdd)

        Mnotinv = N[1 0; 2 0]
        @assert !isinvertible(Mnotinv)

        Minv = N[1 2; -1 0.4]
        @assert isinvertible(Minv)

        # invertible matrix times a bounded polyhedron
        L = linear_map(Minv, Pbdd)
        @test L isa HPolyhedron{N}

        # invertible matrix times an unbounded polyhedron
        L = linear_map(Minv, Punbdd)
        @test L isa HPolyhedron{N}

        # not invertible matrix times a bounded polyhedron
        L = linear_map(Mnotinv, Pbdd, algorithm="vrep") # Requires Polyhedra because it works on vertices
        @test L isa VPolytope

        L = linear_map(Mnotinv, Pbdd, algorithm="elimination")
        @test L isa HPolyhedron

        # test default
        L = linear_map(Mnotinv, Pbdd)
        @test L isa HPolyhedron

        # not invertible matrix times an unbounded polyhedron
        @test linear_map(Mnotinv, Punbdd) isa HPolyhedron

        # check that we can use sparse matrices as well ; Requires SparseArrays
        L = linear_map(sparse(Minv), Pbdd, algorithm="inv_right")
        @test L isa HPolyhedron{N}
        L = linear_map(sparse(Minv), Punbdd, algorithm="inv_right")
        @test L isa HPolyhedron{N}
        L = linear_map(sparse(Mnotinv), Pbdd, algorithm="vrep") # Requires Polyhedra because it works on vertices
        @test L isa VPolytope
        # breaks because "inv_right" requires an invertible matrix
        @test_throws ArgumentError linear_map(sparse(Mnotinv), Punbdd, algorithm="inv_right")

        # remove a repeated constraint (#909)
        Q = HPolyhedron([HalfSpace([-1.29817, 1.04012], 6.07731),
                         HalfSpace([-1.29348, -0.0920708], 1.89515)])
        addconstraint!(Q, Q.constraints[2])
        remove_redundant_constraints!(Q)
        @test length(constraints_list(Q)) == 2

        # concrete projection of an unbounded set
        # P = {x, y, z : x >= 0, y >= 0, z >= 0} is the positive orthant
        P = HPolyhedron([HalfSpace(N[-1, 0, 0], N(0)),
                         HalfSpace(N[0, -1, 0], N(0)),
                         HalfSpace(N[0, 0, -1], N(0))])
        πP = project(P, [1, 2])
        @test πP isa HPolyhedron{N}
        @test ispermutation(constraints_list(πP), [HalfSpace(N[-1, 0], N(0)), HalfSpace(N[0, -1], N(0))])

        # projection in unconstrained dimensions
        P = HPolyhedron([HalfSpace(N[1, 0], N(1)), HalfSpace(N[-1, 0], N(0))])
        πP = project(P, [2])
        @test πP isa Universe{N} && dim(πP) == 2
    end

    # tests that require Symbolics
    @static if isdefined(@__MODULE__, :Symbolics)
        vars = @variables x y
        p1 = HPolyhedron([x + y <= 1, x + y >= -1,  x - y <= 1, x - y >= -1], vars)
        b1 = Ball1(zeros(2), 1.0)
        @test isequivalent(p1, b1)

        p2 = HPolyhedron([x == 0, y <= 0], vars)
        h2 = HPolyhedron([HalfSpace([1.0, 0.0], 0.0), HalfSpace([-1.0, 0.0], 0.0), HalfSpace([0.0, 1.0], 0.0)])
        @test p2 ⊆ h2 && h2 ⊆ p2 # isequivalent(p2, h2) see #2370
    end
end
