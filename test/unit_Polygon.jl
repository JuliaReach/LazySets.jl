for N in [Float64, Float32, Rational{Int}]
    # Empty polygon
    p = HPolygon{N}()

    # Test that constraints are automatically sorted
    c1 = LinearConstraint(N[2., 2.], N(12.))
    c2 = LinearConstraint(N[-3., 3.], N(6.))
    c3 = LinearConstraint(N[-1., -1.], N(0.))
    c4 = LinearConstraint(N[2., -4.], N(0.))
    addconstraint!(p, c3)
    addconstraint!(p, c1)
    addconstraint!(p, c4)
    addconstraint!(p, c2)
    @test p.constraints[1] == c1
    @test p.constraints[2] == c2
    @test p.constraints[3] == c3
    @test p.constraints[4] == c4

    # conversion to optimized polygon
    po = HPolygonOpt(p)
    # conversion back
    HPolygon(po)
    # conversion from HPolytope
    polytope = HPolytope{N}()
    addconstraint!(polytope, c1)
    addconstraint!(polytope, c2)
    addconstraint!(polytope, c3)
    addconstraint!(polytope, c4)
    HPolygon(polytope)
    HPolygonOpt(polytope)
    # conversion to HPolytope
    HPolytope(p)
    HPolytope(po)

    # support vector of empty polygon
    @test_throws ErrorException σ([0.], HPolygon())
    @test_throws ErrorException σ([0.], HPolygonOpt(HPolygon()))
    @test_throws ErrorException σ([0.], HPolytope())

    # HPolygon/HPolygonOpt tests
    for p in [p, po]
        # Test Dimension
        @test dim(p) == 2

        # Test Support Vector
        d = N[1., 0.]
        @test σ(d, p) == N[4., 2.]
        d = N[0., 1.]
        @test σ(d, p) == N[2., 4.]
        d = N[-1., 0.]
        @test σ(d, p) == N[-1., 1.]
        d = N[0., -1.]
        @test σ(d, p) == N[0., 0.]
        d = N[1., -1.]
        @test σ(d, p) == N[4., 2.]

        # membership
        @test ∈(N[0., 0.], p)
        @test ∈(N[4., 2.], p)
        @test ∈(N[2., 4.], p)
        @test ∈(N[-1., 1.], p)
        @test ∈(N[2., 3.], p)
        @test ∈(N[1., 1.], p)
        @test ∈(N[3., 2.], p)
        @test ∈(N[5./4., 7./4.], p)
        @test !∈(N[4., 1.], p)
        @test !∈(N[5., 2.], p)
        @test !∈(N[3., 4.], p)
        @test !∈(N[-1., 5.], p)
        @test !∈(N[-1., 0.], p, N(0.6))
        @test ∈(N[-1., 0.], p, N(0.8))

        # an_element function
        @test an_element(p) ∈ p
        p_shallow = HPolygon{N}()
        @test_throws ErrorException an_element(p_shallow)
        addconstraint!(p_shallow, c1)
        @test_throws ErrorException an_element(p_shallow)

        # hrep conversion
        @test tohrep(p) == p

        # constraints list getter
        @test constraints_list(p) == p.constraints
    end

    # Test VRepresentation
    vp = tovrep(p)
    @test N[2., 4.] ∈ vp.vertices_list
    @test N[-1., 1.] ∈ vp.vertices_list
    @test N[0., 0.] ∈ vp.vertices_list
    @test N[4., 2.] ∈ vp.vertices_list

    # test convex hull of a set of points using the default algorithm
    points = to_N(N, [[0.9,0.2], [0.4,0.6], [0.2,0.1], [0.1,0.3], [0.3,0.28]])
    vp = VPolygon(points) # by default, a convex hull is run
    @test vertices_list(vp) == to_N(N, [[0.1,0.3], [0.2,0.1], [0.9,0.2], [0.4,0.6] ])

    vp = VPolygon(points, apply_convex_hull=false) # we can turn it off
    @test vertices_list(vp) == to_N(N, [[0.9,0.2], [0.4,0.6], [0.2,0.1], [0.1,0.3], [0.3,0.28]])

    # test for pre-sorted points
    points = to_N(N, [[0.1, 0.3], [0.2, 0.1], [0.4, 0.3], [0.4, 0.6], [0.9, 0.2]])
    vp = VPolygon(points, algorithm="monotone_chain_sorted")
    @test vertices_list(vp) == to_N(N, [[0.1, 0.3], [0.2, 0.1], [0.9, 0.2], [0.4, 0.6]])

    # test support vector of a VPolygon
    d = N[1., 0.]
    @test σ(d, vp) == points[5]
    d = N[0., 1.]
    @test σ(d, vp) == points[4]
    d = N[-1., 0.]
    @test σ(d, vp) == points[1]
    d = N[0., -1.]
    @test σ(d, vp) == points[2]

    # test that #83 is fixed
    v = VPolygon(to_N(N, [[2.0, 3.0]]))
    @test N[2.0, 3.0] ∈ v
    @test !(N[3.0, 2.0] ∈ v)
    v = VPolygon(to_N(N, [[2.0, 3.0], [-1.0, -3.4]]))
    @test to_N(N, [-1.0, -3.4]) ∈ v

    # an_element function
    v = VPolygon(to_N(N, [[2., 3.]]))
    @test an_element(v) ∈ v

    # subset
    p1 = VPolygon(to_N(N, [[0., 0.], [2., 0.]]))
    p2 = VPolygon(to_N(N, [[1., 0.]]))
    @test ⊆(p2, p1) && ⊆(p2, p1, true)[1]
    @test ⊆(HPolygon{N}(), p1)
    @test ⊆(p1, BallInf(N[1., 0.], N(1.)))

    v1 = N[1.0, 0.0]
    v2 = N[1.0, 1.0]
    v3 = N[0.0, 1.0]
    v4 = N[-1.0, 1.0]
    v5 = N[-1.0, 0.0]
    v6 = N[-1.0, -1.0]
    v7 = N[0, -1.0]
    v8 = N[1.0, -1.0]
    v = [v1, v2, v3, v4, v5, v6, v7, v8]

    for (i, vi) in enumerate(v)
        for j in i:8
                @test vi <= v[j]
        end
    end
end

# default Float64 constructor
HPolygon()
