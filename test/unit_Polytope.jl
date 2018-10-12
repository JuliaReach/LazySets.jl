global test_suite_polyhedra

for N in [Float64, Rational{Int}, Float32]
    # -----
    # H-rep
    # -----

    # constructor from matrix and vector
    A = [N(1) N(2); N(-1) N(1)]
    b = [N(1), N(2)]
    p = HPolytope(A, b)
    c = p.constraints
    @test c isa Vector{LinearConstraint{N}}
    @test c[1].a == N[1, 2] && c[1].b == N(1)
    @test c[2].a == N[-1, 1] && c[2].b == N(2)

    # convert back to matrix and vector
    A2, b2 = tosimplehrep(p)
    @test A == A2 && b == b2

    # 2D polytope
    p = HPolytope{N}()
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
    d = N[0, 1]
    @test σ(d, p) == N[2, 4]
    d = N[-1, 0]
    @test σ(d, p) == N[-1, 1]
    d = N[0, -1]
    @test σ(d, p) == N[0, 0]

    # membership
    @test ∈(N[5 / 4, 7 / 4], p)
    @test !∈(N[4, 1], p)

    # singleton list (only available with Polyhedra library)
    if test_suite_polyhedra
        @test length(singleton_list(p)) == 4
    else
        @test_throws AssertionError singleton_list(p)
    end

    if test_suite_polyhedra
        # conversion to and from Polyhedra's VRep data structure
        cl = constraints_list(HPolytope(polyhedron(p)))
        @test length(p.constraints) == length(cl)
    end

    if test_suite_polyhedra
        # convert hyperrectangle to a HPolytope
        H = Hyperrectangle(N[1, 1], N[2, 2])
        P = convert(HPolytope, H)
        vlist = vertices_list(P)
        @test ispermutation(vlist, [N[3, 3], N[3, -1], N[-1, -1], N[-1, 3]])
    end

    # -----
    # V-rep
    # -----

    # constructor from a VPolygon
    polygon = VPolygon([N[0, 0], N[1, 0], N[0, 1]])
    p = VPolytope(polygon)
    @test vertices_list(polygon) == vertices_list(p)

    # dim
    @test dim(p) == 2

    # vertices_list function
    @test vertices_list(p) == p.vertices

    if test_suite_polyhedra
        V = VPolytope(polyhedron(p))

        # conversion to and from Polyhedra's VRep data structure
        vl = vertices_list(V)
        @test length(p.vertices) == length(vl) && vl ⊆ p.vertices

        # list of constraints of a VPolytope; calculates tohrep
        @test ispermutation(constraints_list(V), constraints_list(tohrep(V)))
    end
end

# default Float64 constructors
empty_H_polytope = HPolytope()
empty_V_polytope = VPolytope()
@test empty_H_polytope isa HPolytope{Float64}
@test empty_V_polytope isa VPolytope{Float64}

if test_suite_polyhedra
    # check isempty function of the AbstractPolytope interface
    @test !isempty(empty_H_polytope)
    @test isempty(empty_V_polytope)
end

# Polyhedra tests that only work with Float64
if test_suite_polyhedra
    for N in [Float64]
        # -----
        # H-rep
        # -----

        # intersection
        A = [N(0) N(1); N(1) N(0); N(2) N(2)]
        b = N[0, 0, 1]
        p1 = HPolytope(A, b)
        A = [N(0) N(-1); N(-1) N(0); N(1) N(1)]
        b = N[-0.25, -0.25, -0]
        p2 = HPolytope(A, b)
        cap = intersection(p1, p2)
        # currently broken, see #565

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
        p1 = HPolytope(A, b)
        p2 = HPolytope(A, b)
        cp = cartesian_product(p1, p2)
        cl = constraints_list(cp)
        @test length(cl) == 4

        # vertices_list
        A = [N(1) N(-1)]'
        b = N[1, 0]
        p = HPolytope(A, b)
        vl = vertices_list(p)
        @test ispermutation(vl, [N[0], N[1]])

        # tovrep from HPolytope
        A = [N(0) N(-1); N(-1) N(0); N(1) N(1)]
        b = N[-0.25, -0.25, -0]
        P = HPolytope(A, b)
        @test tovrep(P) isa VPolytope
        @test tohrep(P) isa HPolytope # test no-op

        # checking for emptiness
        P = HPolytope([LinearConstraint(N[1, 0], N(0))])    # x <= 0
        @test !isempty(P)
        addconstraint!(P, LinearConstraint(N[-1, 0], N(-1)))  # x >= 1
        @test isempty(P)

        # -----
        # V-rep
        # -----

        # support vector (only available with Polyhedra library)
        polygon = VPolygon([N[0, 0], N[1, 0], N[0, 1]])
        p = VPolytope(polygon)
        d = N[1, 0]
        @test_throws ErrorException σ(d, p, algorithm="xyz")
        if test_suite_polyhedra
            @test σ(d, p) == N[1, 0]
            p_unbounded = HPolytope([LinearConstraint(N[-1, 0], N(0))])
            @test_throws ErrorException σ(d, p_unbounded)
            p_infeasible = HPolytope([LinearConstraint(N[1], N(0)),
                                      LinearConstraint(N[-1], N(-1))])
            @test_throws ErrorException σ(N[1], p_infeasible)
        else
            @test_throws AssertionError σ(d, p)
        end

        # intersection
        p1 = VPolytope(vertices_list(BallInf(N[0, 0], N(1))))
        p2 = VPolytope(vertices_list(BallInf(N[2, 2], N(1))))
        cap = intersection(p1, p2)
        @test vertices_list(cap) == [N[1, 1]]
        # other polytopic sets
        p3 = VPolygon(vertices_list(p2))
        cap = intersection(p1, p3)
        @test vertices_list(cap) == [N[1, 1]]
        p4 = BallInf(N[2, 2], N(1))
        cap = intersection(p1, p4)
        @test vertices_list(cap) == [N[1, 1]]

        # convex hull
        v1 = N[1, 0]
        v2 = N[0, 1]
        v3 = N[-1, 0]
        v4 = N[0, -1]
        v5 = N[0, 0]
        p1 = VPolytope([v1, v2, v5])
        p2 = VPolytope([v3, v4])
        ch = convex_hull(p1, p2)
        vl = vertices_list(ch)
        @test v1 ∈ vl && v2 ∈ vl && v3 ∈ vl && v4 ∈ vl
        # Note: The redundant vertex v5 is not removed (see #561).
        # This test can be removed (and the length above should be corrected)
        # when that issue is resolved.
        @test length(vl) == 5 && v5 ∈ vl

        # Cartesian product
        p1 = VPolytope([N[0, 0], N[1, 1]])
        p2 = VPolytope([N[2]])
        cp = cartesian_product(p1, p2)
        vl = vertices_list(cp)
        @test ispermutation(vl, [N[0, 0, 2], N[1, 1, 2]])

        # tohrep from VPolytope
        P = VPolytope([v1, v2, v3, v4, v5])
        @test tohrep(P) isa HPolytope
        @test tovrep(P) isa VPolytope # no-op
    end
end
