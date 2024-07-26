using LazySets: ⪯

for N in [Float64, Float32, Rational{Int}]
    # Empty polygon
    p = HPolygon{N}()

    # Test that constraints are automatically sorted
    c1 = LinearConstraint(N[2, 2], N(12))
    c2 = LinearConstraint(N[-3, 3], N(6))
    c3 = LinearConstraint(N[-1, -1], N(0))
    c4 = LinearConstraint(N[2, -4], N(0))
    addconstraint!(p, c3)
    addconstraint!(p, c1)
    addconstraint!(p, c4)
    addconstraint!(p, c2)
    @test p.constraints[1] == c1
    @test p.constraints[2] == c2
    @test p.constraints[3] == c3
    @test p.constraints[4] == c4

    # redundancy of constraints
    h1 = HalfSpace([N(1), N(1)], N(1))
    h2 = HalfSpace([N(0), N(1)], N(1)) # redundant
    h3 = HalfSpace([N(-1), N(1)], N(1))
    h4 = HalfSpace([N(0), N(-1)], N(0))
    @test LazySets.isredundant(h2, h1, h3)
    @test !LazySets.isredundant(h3, h2, h4)
    @test LazySets.isredundant(h2, h2, h3) && LazySets.isredundant(h3, h2, h3)
    p2 = HPolygon([h1, h2, h3, h4]; sort_constraints=false) # sorted in the right order
    @test length(p2.constraints) == 4
    remove_redundant_constraints!(p2)
    @test length(p2.constraints) == 3
    h1 = HalfSpace([N(1), N(0)], N(1))
    h2 = HalfSpace([N(0), N(1)], N(1))
    h3 = HalfSpace([N(-1), N(0)], N(1))
    h4 = HalfSpace([N(0), N(-1)], N(1))
    p2 = HPolygon{N}()
    addconstraint!(p2, h1)
    addconstraint!(p2, h2)
    addconstraint!(p2, h3)
    addconstraint!(p2, h4)
    addconstraint!(p2, h1)
    @test ispermutation(p2.constraints, [h1, h2, h3, h4])
    h5 = HalfSpace([N(1), N(1)], N(0))
    addconstraint!(p2, h5)
    @test ispermutation(p2.constraints, [h3, h4, h5])

    # constructor from simple H-rep
    A = N[2 0; 1 3]
    b = N[-1, 1]
    @test constraints_list(HPolygon(A, b)) ==
          constraints_list(HPolygonOpt(A, b)) ==
          [HalfSpace(N[2, 0], N(-1)), HalfSpace(N[1, 3], N(1))]

    # conversion to optimized polygon
    po = convert(HPolygonOpt, p)
    # conversion back
    convert(HPolygon, po)
    # conversion from HPolytope
    polytope = HPolytope{N}()
    addconstraint!(polytope, c1)
    addconstraint!(polytope, c2)
    addconstraint!(polytope, c3)
    addconstraint!(polytope, c4)
    convert(HPolygon, polytope)
    convert(HPolygonOpt, polytope)
    # conversion to HPolytope
    HPolytope(constraints_list(p))
    HPolytope(constraints_list(po))

    # conversion from other set type
    H = Hyperrectangle(; low=N[-1, -1], high=N[1, 1])
    HPolygon(constraints_list(H))
    HPolygonOpt(constraints_list(H))

    # support vector of polygon with no constraints
    @test_throws AssertionError σ(N[0], HPolygon{N}())
    @test_throws AssertionError σ(N[0], HPolygonOpt{N}())

    # boundedness
    @test isbounded(p) && isbounded(po)
    @test !isbounded(HPolygon{N}(), false) &&
          !isbounded(HPolygonOpt{N}(), false)
    @test_throws AssertionError HPolygon(LinearConstraint{N,Vector{N}}[];
                                         check_boundedness=true)
    @test_throws AssertionError HPolygonOpt(LinearConstraint{N,Vector{N}}[];
                                            check_boundedness=true)

    # subset
    b1 = Ball1(N[0, 0], N(1))
    b2 = Ball1(N[1, 1], N(4))
    l1 = LinearMap(N[1 0; 0 1], b1)
    l2 = LinearMap(N[1 0; 0 1], b2)
    subset, point = ⊆(p, b1, true)
    @test !subset && point ∈ p && point ∉ b1
    subset, point = ⊆(p, l1, true)
    @test !subset && point ∈ p && point ∉ l1
    subset, point = ⊆(p, b2, true)
    @test subset && p ⊆ b2 && point == N[]
    subset, point = ⊆(p, l2, true)
    @test subset && p ⊆ l2 && point == N[]
    @test p ⊆ p

    # HPolygon/HPolygonOpt tests
    for (hp, t_hp) in [(p, HPolygon), (po, HPolygonOpt)]
        # constructors
        @test t_hp{N}() isa t_hp{N}
        @test t_hp{N,Vector{N}}() isa t_hp{N,Vector{N}}
        clist = [HalfSpace(N[1, 0], N(1)), HalfSpace(sparsevec([1], N[-1], 2), N(-1))]
        P = t_hp(clist)
        @test P isa t_hp{N,Vector{N}} &&
              P.constraints == [HalfSpace(N[1, 0], N(1)), HalfSpace(N[-1, 0], N(-1))]

        # Test Dimension
        @test dim(hp) == 2

        # test support vector, with linear and binary search
        d = N[1, 0]
        @test σ(d, hp) == N[4, 2]
        @test σ(d, hp; linear_search=true) == σ(d, hp; linear_search=false)
        d = N[0, 1]
        @test σ(d, hp) == N[2, 4]
        @test σ(d, hp; linear_search=true) == σ(d, hp; linear_search=false)
        d = N[0, -1]
        @test σ(d, hp) == N[0, 0]
        @test σ(d, hp; linear_search=true) == σ(d, hp; linear_search=false)
        d = N[-1, 0]
        @test σ(d, hp) == N[-1, 1]
        @test σ(d, hp; linear_search=true) == σ(d, hp; linear_search=false)
        d = N[1, -1]
        @test σ(d, hp) == N[4, 2]
        @test σ(d, hp; linear_search=true) == σ(d, hp; linear_search=false)

        # membership
        @test N[0, 0] ∈ hp
        @test N[4, 2] ∈ hp
        @test N[2, 4] ∈ hp
        @test N[-1, 1] ∈ hp
        @test N[2, 3] ∈ hp
        @test N[1, 1] ∈ hp
        @test N[3, 2] ∈ hp
        @test N[5 / 4, 7 / 4] ∈ hp
        @test N[4, 1] ∉ hp
        @test N[5, 2] ∉ hp
        @test N[3, 4] ∉ hp
        @test N[-1, 5] ∉ hp

        # an_element function
        @test an_element(hp) ∈ hp
        hp_shallow = HPolygon{N}()
        @test_throws AssertionError an_element(hp_shallow)
        addconstraint!(hp_shallow, c1)
        @test_throws AssertionError an_element(hp_shallow)

        # ispolyhedral
        @test ispolyhedral(hp)

        # isuniversal
        answer, w = isuniversal(hp, true)
        @test !isuniversal(hp) && !answer && w ∉ hp

        # isempty
        @test !isempty(hp)

        # hrep conversion
        @test tohrep(hp) == hp

        # constraints list getter
        @test constraints_list(hp) == hp.constraints

        # test constraints list of a VPolygon
        vp = tovrep(hp)
        @test ispermutation(constraints_list(vp), hp.constraints)

        # translation
        @test translate(hp, N[1, 2]) == t_hp([HalfSpace(N[2, 2], N(18)), HalfSpace(N[-3, 3], N(9)),
                                              HalfSpace(N[-1, -1], N(-3)), HalfSpace(N[2, -4], N(-6))])

        # test for concrete minkowski sum
        A = [N[4, 0], N[6, 2], N[4, 4]]
        B = [N[-2, -2], N[2, 0], N[2, 2], N[-2, 4]]

        P = VPolygon(A)
        Q = VPolygon(B)
        PQ = minkowski_sum(P, Q)
        @test is_cyclic_permutation(PQ.vertices,
                                    [N[2, -2], N[6, 0], N[8, 2],
                                     N[8, 4], N[6, 6], N[2, 8]])

        # test for different starting points in vertices_list of minkowski sum
        P = VPolygon([N[4, 0], N[6, 2], N[4, 4]])
        P2 = VPolygon([N[4, 4], N[4, 0], N[6, 2]])
        Q = VPolygon([N[-2, -2], N[2, 0], N[2, 2], N[-2, 4]])
        @test is_cyclic_permutation(vertices_list(minkowski_sum(P, Q)),
                                    vertices_list(minkowski_sum(P2, Q)))

        # test for corner case of parallel edges in minkowski sum
        C = [N[10, 5], N[10, 10], N[8, 10], N[5, 8], N[5, 5]]
        D = [N[-1, -2], N[1, -2], N[-1, 2]]

        R = VPolygon(C)
        S = VPolygon(D)
        RS = minkowski_sum(R, S)
        @test is_cyclic_permutation(RS.vertices,
                                    [N[4, 3], N[11, 3], N[11, 8],
                                     N[9, 12], N[7, 12], N[4, 10]])
    end

    # low/high/extrema
    p2 = VPolygon([N[5, 1], N[4, 0], N[3, 1], N[4, 2]])
    @test low(p2) == N[3, 0] && low(p2, 1) == N(3) && low(p2, 2) == N(0)
    @test high(p2) == N[5, 2] && high(p2, 1) == N(5) && high(p2, 2) == N(2)
    @test extrema(p2) == (N[3, 0], N[5, 2]) && extrema(p2, 1) == (N(3), N(5)) &&
          extrema(p2, 2) == (N(0), N(2))
    # singleton
    p2 = VPolygon([N[1, 2]])
    @test low(p2) == N[1, 2] && low(p2, 1) == N(1) && low(p2, 2) == N(2)
    @test high(p2) == N[1, 2] && high(p2, 1) == N(1) && high(p2, 2) == N(2)
    @test extrema(p2) == (N[1, 2], N[1, 2]) && extrema(p2, 1) == (N(1), N(1)) &&
          extrema(p2, 2) == (N(2), N(2))
    # empty polygon
    p2 = VPolygon{N}()
    @test low(p2) == N[Inf, Inf] && low(p2, 1) == N(Inf) && low(p2, 2) == N(Inf)
    @test high(p2) == N[-Inf, -Inf] && high(p2, 1) == N(-Inf) && high(p2, 2) == N(-Inf)
    @test extrema(p2) == (N[Inf, Inf], N[-Inf, -Inf]) && extrema(p2, 1) == (N(Inf), N(-Inf)) &&
          extrema(p2, 2) == (N(Inf), N(-Inf))

    # ===================================
    # Concrete intersection
    # ===================================

    # concrete intersection of H-rep
    p2 = HPolygon{N}()
    c1 = LinearConstraint(N[2, 2], N(14))
    c2 = LinearConstraint(N[-3, 3], N(6))
    c3 = LinearConstraint(N[-1, -1], N(-2))
    c4 = LinearConstraint(N[2, -4], N(0))
    addconstraint!(p2, c3)
    addconstraint!(p2, c1)
    addconstraint!(p2, c4)
    addconstraint!(p2, c2)
    p3 = intersection(p, p2)
    @test length(constraints_list(p3)) == 4

    # concrete intersection of V-rep
    paux = VPolygon([N[0, 0], N[1, 0], N[0, 1], N[1, 1]])
    qaux = VPolygon([N[1, -1 / 2], N[-1 / 2, 1], N[-1 / 2, -1 / 2]])
    xaux = intersection(paux, qaux)
    oaux = VPolygon([N[0, 0], N[1 / 2, 0], N[0, 1 / 2]])
    @test xaux ⊆ oaux && oaux ⊆ xaux # TODO use isequivalent
    @test LazySets._intersection_vrep_2d(paux.vertices, qaux.vertices) == xaux.vertices
    # empty intersection
    qaux = VPolygon([N[2, 2]])
    xaux = intersection(paux, qaux)
    @test xaux == EmptySet{N}(2)

    # concrete intersection of H-rep and V-rep
    q1 = intersection(paux, p)
    q2 = intersection(p, paux)
    @test isequivalent(q1, q2)
    # HPolygons have the same constraints in the same order
    @test q1 isa HPolygon && q2 isa HPolygon && q1 == q2

    # check that tighter constraints are used in intersection (#883)
    h1 = HalfSpace([N(1), N(0)], N(3))
    h2 = HalfSpace([N(-1), N(0)], N(-3))
    h3 = HalfSpace([N(0), N(1)], N(7))
    h4 = HalfSpace([N(0), N(-1)], N(-5))
    h5 = HalfSpace([N(0), N(1)], N(6))
    h6 = HalfSpace([N(0), N(-1)], N(-4))
    p1 = HPolygon([h1, h3, h2, h4])
    p2 = HPolygon([h1, h5, h2, h6])
    c = intersection(p1, p2).constraints
    @test c == [h1, h5, h2, h4]
    h1 = HalfSpace([N(1), N(1)], N(1))
    h2 = HalfSpace([N(-1), N(1)], N(1))
    h3 = HalfSpace([N(0), N(-1)], N(0))
    p1 = HPolygon([h1, h2, h3])
    h4 = HalfSpace([N(0), N(1)], N(1))
    h5 = HalfSpace([N(-1), N(-1)], N(0))
    h6 = HalfSpace([N(1), N(-1)], N(0))
    p2 = HPolygon([h4, h5, h6])
    c = intersection(p1, p2).constraints
    @test c == [h1, h2, h5, h6]

    # check that empty polygon (infeasible constraints) has no vertices (#918)
    P = HPolygon([HalfSpace(N[1, 1], N(0)), HalfSpace(N[-1, 0], N(-1)),
                  HalfSpace(N[0, -1], N(-1))])
    @test vertices_list(P) == Vector{Vector{N}}()
    P = HPolygon([HalfSpace(N[1, 0], N(0)), HalfSpace(N[-1, 0], N(-1))])
    @test vertices_list(P) == Vector{Vector{N}}()
    @test isempty(P)

    # empty intersection results in empty set
    p3 = tohrep(VPolygon([N[0, 0]]))
    p4 = tohrep(VPolygon([N[1, 1]]))
    @test intersection(p3, p4) isa EmptySet{N}

    # concrete linear map
    # in 2D and for an invertible map we get an HPolygon (see #631 and #1093)
    HP = convert(HPolygon, BallInf(N[0, 0], N(1)))
    @test linear_map(N[1 0; 0 2], HP) isa HPolygon{N}
    # in higher dimensions we get an HPolytope (#2168)
    if test_suite_polyhedra
        p4 = convert(HPolygon, BallInf(zeros(N, 2), N(1)))
        A = ones(N, 4, 2)
        @test linear_map(A, p4) isa HPolytope
    end
    # linear_map with redundant vertices
    A = N[1 0; 0 0]
    P = VPolygon([N[1, 1], N[-1, 1], N[1, -1], N[-1, -1]])
    Q = linear_map(A, P; apply_convex_hull=true)
    @test ispermutation(vertices_list(Q), [N[1, 0], N[-1, 0]])
    # non-square matrices
    A = N[1 0]
    Q = linear_map(A, P)
    @test !(Q isa VPolygon) && dim(Q) == 1
    A = N[1 0; 0 1; 1 1]
    Q = linear_map(A, P)
    @test !(Q isa VPolygon) && dim(Q) == 3

    # vertices_list removes duplicates by default (#1405)
    p3 = HPolygon([HalfSpace(N[1, 0], N(0)), HalfSpace(N[0, 1], N(0)),
                   HalfSpace(N[-1, 0], N(0)), HalfSpace(N[0, -1], N(0))])
    @test vertices_list(p3) == [N[0, 0]]

    # Test VRepresentation
    vp = tovrep(p)
    @test N[2, 4] ∈ vertices_list(vp)
    @test N[-1, 1] ∈ vertices_list(vp)
    @test N[0, 0] ∈ vertices_list(vp)
    @test N[4, 2] ∈ vertices_list(vp)
    @test tovrep(vp) == vp

    # test convex hull of a set of points using the default algorithm
    v1 = N[1 // 10, 3 // 10]
    v2 = N[1 // 5, 1 // 10]
    v3 = N[2 // 5, 3 // 5]
    v4 = N[9 // 10, 1 // 5]
    v5 = N[3 // 10, 7 // 25]
    points = [v4, v3, v2, v1, v5]
    vp = VPolygon(points) # by default, a convex hull is run
    @test ispermutation(vertices_list(vp), [v1, v2, v3, v4])

    vp = VPolygon(points; apply_convex_hull=false) # we can turn it off
    @test ispermutation(vertices_list(vp), [v1, v2, v3, v4, v5])

    # test for pre-sorted points
    v3_new = N[2 // 5, 3 // 10]
    points = [v1, v2, v3_new, v3, v4]
    vp = VPolygon(points; algorithm="monotone_chain_sorted")
    @test ispermutation(vertices_list(vp), [v1, v2, v3, v4])

    # test that #83 is fixed
    p = VPolygon([N[2, 3]])
    @test N[2, 3] ∈ p
    @test N[3, 2] ∉ p
    v = N[-1, -17 // 5]
    p = VPolygon([N[2, 3], v])
    @test v ∈ p

    P = VPolygon([N[0, 0], [0, 1], [0, 2]]; apply_convex_hull=false)
    Q = remove_redundant_vertices(P)
    @test length(P.vertices) == 3 && length(Q.vertices) == 2

    # ispolyhedral
    @test ispolyhedral(p)

    # an_element function
    p = VPolygon([N[2, 3]])
    @test an_element(p) ∈ p

    # membership
    point = N[0, 1]
    @test point ∉ VPolygon(Vector{Vector{N}}(); apply_convex_hull=false)
    @test point ∉ VPolygon([N[0, 2]])
    @test point ∈ VPolygon([N[0, 0], N[0, 2]])
    @test point ∉ VPolygon([N[1, 0], N[1, 2]])
    @test point ∈ VPolygon([N[1, 1], N[0, 0], N[-1, 1]])
    @test point ∉ VPolygon([N[1, 1], N[0, 0], N[1, 0]])

    # translation
    vp2 = VPolygon([N[0, 0], N[1, 0], N[0, 1]])
    @test translate(vp2, N[1, 2]) == VPolygon([N[1, 2], N[2, 2], N[1, 3]])
    vp2cp = copy(vp2)
    @test translate!(vp2cp, N[1, 2]) == VPolygon([N[1, 2], N[2, 2], N[1, 3]]) == vp2cp

    # subset
    p1 = VPolygon([N[0, 0], N[2, 0]])
    p2 = VPolygon([N[1, 0]])
    b = BallInf(N[2, 0], N(1))
    @test p2 ⊆ p1 && ⊆(p2, p1, true)[1]
    subset, witness = ⊆(p1, b, true)
    @test p1 ⊈ b && !subset && witness ∈ p1 && witness ∉ b
    subset, witness = ⊆(p2, b, true)
    @test p2 ⊆ b && subset

    v1 = N[1, 0]
    v2 = N[1, 1]
    v3 = N[0, 1]
    v4 = N[-1, 1]
    v5 = N[-1, 0]
    v6 = N[-1, -1]
    v7 = N[0, -1]
    v8 = N[1, -1]
    v = [v1, v2, v3, v4, v5, v6, v7, v8]

    for (i, vi) in enumerate(v)
        for j in i:8
            @test vi ⪯ v[j]
        end
    end

    # hrep conversion
    v1 = N[9 // 10, 1 // 5]
    v2 = N[2 // 5, 3 // 5]
    v3 = N[1 // 5, 1 // 10]
    v4 = N[1 // 10, 3 // 10]
    v5 = N[7 // 10, 2 // 5]
    v6 = N[4 // 5, 3 // 10]
    v7 = N[3 // 10, 11 // 20]
    v8 = N[1 // 5, 9 // 20]
    points5 = [v1, v2, v3, v4, v5, v6, v7, v8]
    for i in [0, 1, 2, 4, 8]
        points = i == 0 ? Vector{Vector{N}}() : points5[1:i]
        vp = VPolygon(points; apply_convex_hull=i > 0)
        h1 = tohrep(vp)
        @test convert(HPolygon, vp) == h1
        if i == 0
            @test isempty(h1)
            continue
        elseif i == 1
            @test v1 ∈ h1
        elseif i == 2
            c = h1.constraints[1]
            @test c.a ≈ N[2 // 5, 1 // 2] && c.b ≈ N(23 // 50)
            c = h1.constraints[3]
            @test c.a ≈ N[-2 // 5, -1 // 2] && c.b ≈ N(-23 // 50)
        elseif i == 4
            @test length(h1.constraints) == 4
            c = h1.constraints[1]
            @test c.a ≈ N[2 // 5, 1 // 2] && c.b ≈ N(23 // 50)
            c = h1.constraints[2]
            @test c.a ≈ N[-3 // 10, 3 // 10] && c.b ≈ N(3 // 50)
            c = h1.constraints[3]
            @test c.a ≈ N[-1 // 5, -1 // 10] && c.b ≈ N(-1 // 20)
            c = h1.constraints[4]
            @test c.a ≈ N[1 // 10, -7 // 10] && c.b ≈ N(-1 // 20)
        end

        # check that constraints are sorted correctly
        h2 = HPolygon{N}()
        for c in h1.constraints
            addconstraint!(h2, c)
        end
        @test h1.constraints == h2.constraints
    end

    # empty VPolygon: conversion to hrep
    @test tohrep(VPolygon{N}()) isa EmptySet{N}

    # test VPolygon constructor given the matrix of vertices
    m = N[4 0; 6 2; 4 4]'
    P = VPolygon(m)
    @test is_cyclic_permutation(vertices_list(P), [N[4, 0], N[6, 2], N[4, 4]])

    # test concrete projection
    V = VPolygon([N[0, 1], N[1, 0], N[-1, 0]])
    @test project(V, [1]) == Interval(N(-1), N(1))
    @test project(V, 1:2) == V
    @test project(V, [2, 1]) == VPolygon([N[1, 0], N[0, 1], N[0, -1]])
    V = VPolygon([N[1, 0], N[1, 1]])
    @test project(V, [1]) == Interval(N(1), N(1))
    @test_throws AssertionError project(V, [3])
    @test_throws AssertionError project(V, [1, 3])
    @test_throws ArgumentError project(V, [1, 2, 3])

    # permute
    V = VPolygon([N[1, 0], N[1, 2]])
    @test permute(V, [2, 1]) == VPolygon([N[0, 1], N[2, 1]])

    # concrete cartesian product
    V = VPolygon([N[0, 1], N[1, 0], N[-1, 0]])
    I = Interval(N(0), N(1))
    Vcp = cartesian_product(I, V)
    Vcp′ = VPolytope([N[0, -1, 0], N[0, 1, 0], N[0, 0, 1], N[1, -1, 0], N[1, 1, 0], N[1, 0, 1]])
    @test isequivalent(cartesian_product(I, V), Vcp)
    @test LazySets._issubset_vertices_list(Vcp, Vcp′, false)
    @test LazySets._issubset_vertices_list(Vcp′, Vcp, false)

    # concrete projection of a cartesian product
    @test project(I × V, 2:3) === V
    @test project(I × V, 1:1) == I
    @test isequivalent(project(I × V, 1:2), VPolygon([N[0, -1], N[1, -1], N[1, 1], [0, 1]]))

    # rectification
    P = VPolygon([N[-1, 2], N[1, 1], N[1, 3]])
    rectify(P)
    Q1 = rectify(P, true)
    Q2 = rectify(P, false)
    for (d, res) in [(N[-1, 0], N(0)),
                     (N[0, -1], N(-1)),
                     (N[-1, -1], N(-3 // 2))]
        @test ρ(d, Q1) == res
        @test abs(ρ(d, Q1) - ρ(d, Q2)) < 1e-8  # precision of lazy intersection not good
    end

    # concrete Minkowski sum for arbitrary polytopic sets performs 2D computation
    X = BallInf(zeros(N, 2), N(1))
    Y = Ball1(zeros(N, 2), N(1))
    Z = minkowski_sum(X, Y)
    # the vertices are [±1, ±2] and [±2, ±1]
    vlist = [[[N(-1)^i * k, N(-1)^j * (3 - k)] for i in 1:2, j in 1:2, k in 1:2]...]
    @test Z isa VPolygon{N} && ispermutation(vertices_list(Z), vlist)
    Z = minkowski_sum(X + Y, Singleton(zeros(N, 2)))
    @test Z isa VPolygon{N} && ispermutation(vertices_list(Z), vlist)

    # area & volume
    v1 = N[1, 3]
    v2 = N[2, 1]
    v3 = N[4, 6]
    v4 = N[9, 2]
    P = VPolygon([v1, v2, v3, v4])
    @test area(P) == volume(P) == N(21)
    Q = tohrep(P)
    @test area(Q) == volume(Q) == N(21)

    # Matrix to VPolygon
    M = N[0 1 0; 0 0 1]
    Vs = [N[0, 0], N[1, 0], N[0, 1]]
    P = VPolygon(M)
    @test P == VPolygon(Vs)
    @test eltype(P.vertices) == eltype(Vs)
    # StaticArraysCore.SMatrix to VPolygon
    M = @SMatrix N[0 1 0; 0 0 1]
    Vs = [@SVector[N(0), N(0)], @SVector[N(1), N(0)], @SVector[N(0), N(1)]]
    Q = VPolygon(M)
    @test Q == VPolygon(Vs) == P
    @test eltype(Q.vertices) == eltype(Vs)
end

for N in [Float64, Float32]
    v1 = N[1 // 10, 3 // 10]
    v2 = N[1 // 5, 1 // 10]
    v3 = N[2 // 5, 3 // 10]
    v4 = N[2 // 5, 3 // 5]
    v5 = N[9 // 10, 1 // 5]

    # test support vector of a VPolygon
    points = [v1, v2, v3, v4, v5]
    vp = VPolygon(points; algorithm="monotone_chain_sorted")
    d = N[1, 0]
    @test σ(d, vp) == points[5]
    d = N[0, 1]
    @test σ(d, vp) == points[4]
    d = N[-1, 0]
    @test σ(d, vp) == points[1]
    d = N[0, -1]
    @test σ(d, vp) == points[2]
    dirs = [N[1, 0], N[1, 1 // 2], N[1 // 2, 1 // 2], N[1 // 2, 1], N[0, 1],
            N[-1 // 2, 1], N[-1 // 2, 1 // 2], N[-1, 1 // 2], N[-1, 0], N[1, -1 // 2],
            N[1 // 2, -1 // 2], N[1 // 2, -1], N[0, -1], N[-1 // 2, -1],
            N[-1 // 2, -1 // 2], N[-1, -1 // 2]]
    B = Ball2(zeros(N, 2), N(1))
    P = HPolygon([HalfSpace(di, ρ(di, B)) for di in dirs])
    vlistP = vertices_list(P)
    V = VPolygon(vlistP)
    all(x -> _isapprox(σ(x, V), x), vlistP)

    # test adding constraints, with linear and binary search
    p1 = HPolygon{N}()
    p2 = HPolygon{N}()
    po1 = HPolygonOpt{N}()
    po2 = HPolygonOpt{N}()
    n = 10
    for i in 1:n
        constraint = LinearConstraint(rand(N, 2), rand(N))
        addconstraint!(p1, constraint; linear_search=true)
        addconstraint!(p2, constraint; linear_search=i <= 2)
        addconstraint!(po1, constraint; linear_search=true)
        addconstraint!(po2, constraint; linear_search=i <= 2)
    end
    n = length(p1.constraints)
    all_equal = n == length(p2.constraints) == length(po1.constraints) ==
                length(po2.constraints)
    if all_equal  # may be violated, see #2426
        for i in 1:n
            @test allequal([p1.constraints[i], p2.constraints[i],
                            po1.constraints[i], po2.constraints[i]])
        end
    end

    for (hp, t_hp) in [(p1, HPolygon), (po1, HPolygonOpt)]
        # normalization
        p2 = normalize(hp)
        @test p2 isa t_hp{N}
        for hs in constraints_list(p2)
            @test norm(hs.a) ≈ N(1)
        end
    end

    # check redundancy removal
    p2 = HPolygon([HalfSpace(N[-1.29817, 1.04012], N(6.07731)),
                   HalfSpace(N[-1.29348, -0.0920708], N(1.89515)),
                   HalfSpace(N[-1.42909, -0.347449], N(1.50577)),
                   HalfSpace(N[0.349446, -0.130477], N(-0.575904)),
                   HalfSpace(N[1.50108, -0.339291], N(-2.18331)),
                   HalfSpace(N[2.17022, -0.130831], N(-2.14411))])
    addconstraint!(p2, p2.constraints[2])
    @test length(p2.constraints) == 6

    # correct sorting of constraints in intersection (#2187)
    h = Hyperrectangle(N[0.98069, 0.85020],
                       N[0.00221, 0.00077])
    p = overapproximate(Ball2(N[1, 1], N(15 // 100)), 1e-4)
    @test !isempty(intersection(HPolygon(constraints_list(h)), p))

    # test that concrete minkowski sum of a singleton and a polygon works
    x = VPolygon([N[1, 1]])
    y = VPolygon([N[1.1, 2.2], N[0.9, 2.2], N[0.9, 2.0], N[1.1, 2.0]])
    v = minkowski_sum(x, y)

    a = VPolygon([N[2.1, 3.2], N[1.9, 3.2], N[1.9, 3.0], N[2.1, 3.0]])
    @test v == a
    x = VPolygon([N[1, 1]])
    y = VPolygon([N[2, 1]])
    @test minkowski_sum(x, y) == VPolygon([N[3, 2]])

    # rationalize
    P = VPolygon([N[1, 2], N[3, 4]])
    @test rationalize(P) == VPolygon([[1 // 1, 2 // 1], [3 // 1, 4 // 1]])

    # test corner case with binary search (#2426)
    c1 = HalfSpace(N[0.9636463837316949, 0.8950288819001659], N(0.19114554280256102))
    c2 = HalfSpace(N[0.336341452138065, 0.3367612579000614], N(0.6801677354217885))
    c3 = HalfSpace(N[0.6052463696308072, 0.606001419182304], N(0.4494286571876329))
    for ls in [true, false]
        x = [c1, c2]
        addconstraint!(x, c3; linear_search=ls)
        if N == Float64
            @test x == [c1, c3, c2]
        else
            @test x == [c1, c3]
        end
    end

    # rand
    @test rand(VPolygon; N=N, num_vertices=0) == VPolygon{N}()
    vp = rand(VPolygon; N=N, num_vertices=1)
    @test vp isa VPolygon{N} && length(vp.vertices) == 1
    @test rand(VPolygon; N=N) isa VPolygon{N}
end

for N in [Float64]
    # check boundedness after conversion
    H = Hyperrectangle(; low=N[-1, -1], high=N[1, 1])
    HPolygon(constraints_list(H); check_boundedness=true)
    HPolygonOpt(constraints_list(H); check_boundedness=true)

    p = HPolygon{N}()
    c1 = LinearConstraint(N[2, 2], N(12))
    c2 = LinearConstraint(N[-3, 3], N(6))
    c3 = LinearConstraint(N[-1, -1], N(0))
    c4 = LinearConstraint(N[2, -4], N(0))
    addconstraint!(p, c3)
    addconstraint!(p, c1)
    addconstraint!(p, c4)
    addconstraint!(p, c2)
    po = convert(HPolygonOpt, p)

    # test boundedness
    @test isbounded(p, false) && isbounded(po, false)

    # is intersection empty
    p3 = convert(HPolygon, Hyperrectangle(; low=N[-1, -1], high=N[1, 1]))
    I1 = Interval(N(9 // 10), N(11 // 10))
    I2 = Interval(N(1 // 5), N(3 // 10))
    I3 = Interval(N(4), N(5))
    @test !is_intersection_empty(I1 × I2, p3)
    @test is_intersection_empty(I1 × I3, p3)

    # redundancy with almost-parallel constraints (numeric issues, #2102)
    Hs = [HalfSpace([0.05125918978413134, 0.0], 0.956012965730266),
          HalfSpace([-0.0363879863894044, 0.3036279569180349], -0.370213880061176),
          HalfSpace([0.03529692609786872, -0.2945239519959644], 0.36249837124872736)]
    H = HalfSpace([-0.012707966980287463, 3.809113859067846e-15], -0.2348665397215645)
    addconstraint!(Hs, H)
    @test length(Hs) == 4

    # redundancy with almost-parallel constraints and order issues (#2169)
    Hs = [HalfSpace([1.3877787807814457e-17, 0.0], -7.27918967693285e-18)
          HalfSpace([3.642919299551295e-17, 2.7755575615628914e-17], 3.194076233118193e-17)
          HalfSpace([-6.938893903907228e-17, 0.03196481037863655], 0.1022675780516848)
          HalfSpace([-0.022382796523655775, -3.469446951953614e-18], 0.05268494728581621)
          HalfSpace([-8.673617379884035e-19, -0.01331598581298936], 0.020336411343083463)
          HalfSpace([0.01620615792275367, -2.949029909160572e-17], 0.03845370605895632)]
    P = HPolygon(copy(Hs))
    @test ispermutation(P.constraints, Hs[3:6])

    # redundancy with almost-parallel constraints and approximation issues (#2386)
    P = HPolygon([HalfSpace([1.0, 0.0], 1.0),
                  HalfSpace([0.0, 1.0], 1.0),
                  HalfSpace([-1.0, 0.0], 1.0),
                  HalfSpace([0.0, -1.0], 1.0),
                  HalfSpace([0.17178027783046604, -0.342877222169534], -0.3988040749330524),
                  HalfSpace([0.17228094170254737, -0.3438765582974526], -0.5161575)])
    @test vertices_list(P) == [[-1.0, 1.0]]

    # rationalization
    P = HPolygon([HalfSpace([-1.0, 0.0], 0.0),
                  HalfSpace([0.0, -1.0], 0.0),
                  HalfSpace([1.0, 0.0], 1.0),
                  HalfSpace([0.0, 1.0], 1.0),
                  HalfSpace([-0.2, -0.8], -0.0),
                  HalfSpace([0.6, 0.4], 0.6)])
    @test area(P) ≈ 2 // 3
    @test area(rationalize(P)) == 2 // 3
    Pr = rationalize(BigInt, P, 10 * eps(Float64))
    @test isa(Pr, HPolygon{Rational{BigInt},Vector{Rational{BigInt}}})
end

# random polygons
rand(HPolygon)
rand(HPolygonOpt)

# default Float64 constructors
@test HPolygon() isa HPolygon{Float64,Vector{Float64}}
@test HPolygonOpt() isa HPolygonOpt{Float64,Vector{Float64}}
@test VPolygon() isa VPolygon{Float64}

# isoperationtype
@test !isoperationtype(HPolygon)
@test !isoperationtype(HPolygonOpt)
@test !isoperationtype(VPolygon)
