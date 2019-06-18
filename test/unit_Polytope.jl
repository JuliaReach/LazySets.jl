global test_suite_polyhedra

for N in [Float64, Rational{Int}, Float32]
    # random polytopes
    if test_suite_polyhedra
        rand(HPolytope)
    else
        @test_throws AssertionError rand(HPolytope)
    end
    rand(VPolytope)

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
    @test_throws AssertionError HPolytope(N[1 0; 0 1], N[1, 1];
                                         check_boundedness=true)

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

    # support vector of polytope with no constraints
    @test_throws ErrorException σ(N[0], HPolytope{N}())

    # boundedness
    @test isbounded(p) && isbounded(p, false)
    p2 = HPolytope{N}()
    @test isbounded(p2) && !isbounded(p2, false)

    # membership
    @test N[5 / 4, 7 / 4] ∈ p
    @test N[4, 1] ∉ p

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

    # vertices_list of "universal polytope" (strictly speaking: illegal input)
    @test vertices_list(HPolytope{N}()) == Vector{Vector{N}}()

    if test_suite_polyhedra
        # convert hyperrectangle to a HPolytope
        H = Hyperrectangle(N[1, 1], N[2, 2])
        P = convert(HPolytope, H)
        vlist = vertices_list(P)
        @test ispermutation(vlist, [N[3, 3], N[3, -1], N[-1, -1], N[-1, 3]])
        # check boundedness after conversion
        HPolytope(constraints_list(H); check_boundedness=true)

        # isempty
        @test !isempty(p)
        @test !isempty(HPolytope{N}())

        # H-representaion of an empty v-polytope
        @test tohrep(VPolytope{N}()) == EmptySet{N}()
    end

    # remove redundant constraints
    P = HPolytope([HalfSpace(N[1, 0], N(1)),
                   HalfSpace(N[0, 1], N(1)),
                   HalfSpace(N[-1, -0], N(1)),
                   HalfSpace(N[-0, -1], N(1)),
                   HalfSpace(N[1, 0], N(2))]) # redundant

    Pred = remove_redundant_constraints(P)
    @test length(Pred.constraints) == 4
    @test length(P.constraints) == 5

    # test in-place removal of redundancies
    remove_redundant_constraints!(P)
    @test length(P.constraints) == 4

    # translation
    P2 = translate(P, N[1, 2])
    @test P2 isa HPolytope && ispermutation(constraints_list(P2),
        [HalfSpace(N[1, 0], N(2)), HalfSpace(N[0, 1], N(3)),
         HalfSpace(N[-1, -0], N(0)), HalfSpace(N[-0, -1], N(-1))])

    # subset
    H = BallInf(N[0, 0], N(1))
    P = convert(HPolytope, H)
    @test BallInf(N[0, 0], N(1)) ⊆ P
    @test !(BallInf(N[0, 0], N(1.01)) ⊆ P)

    # =====================
    # Concrete linear map
    # =====================
    LM = linear_map(N[2 3; 1 2], P) # invertible matrix
    @test LM isa HPolytope
    if test_suite_polyhedra
        LM = linear_map(N[2 3; 0 0], P)  # non-invertible matrix
        @test LM isa VPolygon
    end

    # in 4D and for invertible map we get an HPolytope; see #631 and #1093
    HP = convert(HPolytope, Approximations.overapproximate(H * H))
    @test linear_map(Diagonal(N[1, 2, 3, 4]), HP) isa HPolytope

    M = N[2 1; 0 1]
    L1 = linear_map(M, P, use_inv=true)  # calculates inv(M) explicitly
    L2 = linear_map(M, P, use_inv=false) # uses transpose(M) \ c.a for each constraint c of P
    L3 = linear_map(M, P, cond_tol=1e3)  # set a custom tolerance for the condition number (invertibility check)
    # needs M * an_element(Li) to be stable for rational, see #1105
    p = convert(Vector{N}, an_element(P))
    @assert p ∈ P
    @test all([M * p ∈ Li for Li in [L1, L2, L3]])

    # do not check for invertibility => use the vertices
    L4 = linear_map(M, P, check_invertibility=false)
    @test L4 isa VPolygon

    # linear map for mixed types
    M = [2 1; 0 1] # Int's
    LM = linear_map(M, P)
    @test LM isa HPolytope{N}

    # -----
    # V-rep
    # -----

    # constructor from a VPolygon
    polygon = VPolygon([N[0, 0], N[1, 0], N[0, 1]])
    p = VPolytope(polygon)
    @test vertices_list(polygon) == vertices_list(p)

    # dim
    @test dim(p) == 2

    # boundedness
    @test isbounded(p)

    # isempty
    @test !isempty(p)
    @test isempty(VPolytope{N}())

    # vertices_list function
    @test vertices_list(p) == p.vertices

    if test_suite_polyhedra
        V = VPolytope(polyhedron(p))

        # conversion to and from Polyhedra's VRep data structure
        vl = vertices_list(V)
        @test length(p.vertices) == length(vl) && vl ⊆ p.vertices

        # list of constraints of a VPolytope; calculates tohrep
        @test ispermutation(constraints_list(V), constraints_list(tohrep(V)))

        # convert empty VPolytope to a polyhedron
        Vempty = VPolytope()
        @test_throws ErrorException polyhedron(Vempty) # needs to pass the (relative) dim
        Pe = polyhedron(Vempty, relative_dimension=2)
    end

    # membership
    @test N[.49, .49] ∈ p
    @test N[.51, .51] ∉ p

    # translation
    @test translate(p, N[1, 2]) == VPolytope([N[1, 2], N[2, 2], N[1, 3]])

    # copy (see #1002)
    p, q = [N(1)], [N(2)]
    P = VPolytope([p, q])
    Pcopy = copy(P)
    p[1] = N(5)
    # test that Pcopy is independent of P ( = deepcopy)
    @test Pcopy.vertices[1] == [N(1)]
end

# default Float64 constructors
@test HPolytope() isa HPolytope{Float64}
@test VPolytope() isa VPolytope{Float64}

# Polyhedra tests that only work with Float64
if test_suite_polyhedra
    for N in [Float64]
        # -----
        # H-rep
        # -----
        # support function/vector
        d = N[1, 0]
        p_unbounded = HPolytope([LinearConstraint(N[-1, 0], N(0))])
        @test_throws ErrorException σ(d, p_unbounded)
        @test_throws ErrorException ρ(d, p_unbounded)
        p_infeasible = HPolytope([LinearConstraint(N[1], N(0)),
                                  LinearConstraint(N[-1], N(-1))])
        @test_throws ErrorException σ(N[1], p_infeasible)

        # empty intersection
        A = [N(0) N(-1); N(-1) N(0); N(1) N(1)]
        b = N[0, 0, 1]
        p1 = HPolytope(A, b)
        A = [N(0) N(1); N(1) N(0); N(-1) N(-1)]
        b = N[2, 2, -2]
        p2 = HPolytope(A, b)
        @test intersection(p1, p2) isa EmptySet{N}
        @test intersection(p1, p2; backend=Polyhedra.default_library(2, N)) isa EmptySet{N}
        # @test intersection(p1, p2; backend=CDDLib.Library()) isa EmptySet{N}
        # commented because we do not load CDDLib at the moment

        # intersection with half-space
        hs = HalfSpace(N[2, -2], N(-1))
        c1 = constraints_list(intersection(p1, hs))
        c2 = constraints_list(intersection(hs, p1))
        @test length(c1) == 3 && ispermutation(c1, c2)

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

        # checking for emptiness
        P = HPolytope([LinearConstraint(N[1, 0], N(0))])    # x <= 0
        @test !isempty(P)
        addconstraint!(P, LinearConstraint(N[-1, 0], N(-1)))  # x >= 1
        @test isempty(P)

        # checking for empty intersection (also test symmetric calls)
        P = convert(HPolytope, BallInf(zeros(N, 2), N(1)))
        Q = convert(HPolytope, BallInf(ones(N, 2), N(1)))
        R = convert(HPolytope, BallInf(3*ones(N, 2), N(1)))
        res, w = isdisjoint(P, Q, true)
        @test !isdisjoint(P, Q) && !res && w ∈ P && w ∈ Q
        @test !isdisjoint(P, Q; algorithm="sufficient")
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
        res, w = isdisjoint(P, R, true; algorithm="sufficient")
        @test isdisjoint(P, R; algorithm="sufficient") && res && w == N[]

        # test that one can pass a sparse vector as the direction (see #1011)
        P = HPolytope([HalfSpace(N[1, 0], N(1)),
                       HalfSpace(N[0, 1], N(1)),
                       HalfSpace(N[-1, -1], N(-1))])
        @test an_element(P) ∈ P

        # tovrep from HPolytope
        A = N[0 1; 1 0; -1 -1]
        b = N[0.25, 0.25, -0]
        P = HPolytope(A, b)
        @test tovrep(P) isa VPolytope{N}
        @test tohrep(P) isa HPolytope{N}  # no-op

        # -----
        # V-rep
        # -----

        # support vector (only available with Polyhedra library)
        polygon = VPolygon([N[0, 0], N[1, 0], N[0, 1]])
        p = VPolytope(polygon)
        d = N[1, 0]
        @test_throws ErrorException σ(d, p, algorithm="xyz")
        @test σ(d, p) == N[1, 0]

        # intersection
        p1 = VPolytope(vertices_list(BallInf(N[0, 0], N(1))))
        p2 = VPolytope(vertices_list(BallInf(N[2, 2], N(1))))
        cap = intersection(p1, p2)
        @test vertices_list(cap) ≈ [N[1, 1]]
        # other polytopic sets
        p3 = VPolygon(vertices_list(p2))
        cap = intersection(p1, p3)
        @test vertices_list(cap) ≈ [N[1, 1]]
        p4 = BallInf(N[2, 2], N(1))
        cap = intersection(p1, p4)
        vlist = vertices_list(cap)  # contains duplicates (see #1405)
        @test all(v -> v == N[1, 1], vlist)

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
        @test ispermutation(vl, [v1, v2, v3, v4])

        # Cartesian product
        p1 = VPolytope([N[0, 0], N[1, 1]])
        p2 = VPolytope([N[2]])
        cp = cartesian_product(p1, p2)
        vl = vertices_list(cp)
        @test ispermutation(vl, [N[0, 0, 2], N[1, 1, 2]])

        # tohrep from VPolytope
        P = VPolytope([v1, v2, v3, v4, v5])
        @test tohrep(P) isa HPolytope{N}
        @test tovrep(P) isa VPolytope{N}  # no-op

        # subset (see #974)
        H = BallInf(N[0, 0], N(1))
        P = convert(VPolytope, H)
        @test BallInf(N[0, 0], N(1)) ⊆ P
        @test !(BallInf(N[0, 0], N(1.01)) ⊆ P)

        # -----------------
        # mixed H-rep/V-rep
        # -----------------

        # intersection
        p1 = convert(HPolytope, BallInf(N[0, 0], N(1)))
        p2 = convert(VPolytope, BallInf(N[1, 1], N(1)))
        cap1 = intersection(p1, p2)
        cap2 = intersection(p2, p1)
        @test ispermutation([round.(v) for v in vertices_list(cap1)],
                            [N[1, 1], N[0, 1], N[1, 0], N[0, 0]])
        @test ispermutation([round.(v) for v in vertices_list(cap1)],
                            [round.(v) for v in vertices_list(cap2)])
    end
end
