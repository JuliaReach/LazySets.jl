for N in [Float64, Float32] # TODO Rational{Int}
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
    @test p.constraints_list[1] == c1
    @test p.constraints_list[2] == c2
    @test p.constraints_list[3] == c3
    @test p.constraints_list[4] == c4

    # Optimized polygon
    po = HPolygonOpt(p)

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

        # Test containment
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

        # an_element function
        @test an_element(p) ∈ p
    end

    # Test VRepresentation
    vp = tovrep(p)
    @test N[2., 4.] ∈ vp.vertices_list
    @test N[-1., 1.] ∈ vp.vertices_list
    @test N[0., 0.] ∈ vp.vertices_list
    @test N[4., 2.] ∈ vp.vertices_list

    # test convex hull of a set of points using the default algorithm
    points = [N[0.9,0.2], N[0.4,0.6], N[0.2,0.1], N[0.1,0.3], N[0.3,0.28]]
    vp = VPolygon(points) # by default, a convex hull is run
    @test vertices_list(vp) == [ N[0.1,0.3],N[0.2,0.1], N[0.9,0.2],N[0.4,0.6] ]

    vp = VPolygon(points, apply_convex_hull=false) # we can turn it off
    @test vertices_list(vp) == [N[0.9,0.2], N[0.4,0.6], N[0.2,0.1], N[0.1,0.3], N[0.3,0.28]]

    # test for pre-sorted points
    points = [N[0.1, 0.3], N[0.2, 0.1], N[0.4, 0.3], N[0.4, 0.6], N[0.9, 0.2]]
    vp = VPolygon(points, algorithm="monotone_chain_sorted")
    @test vertices_list(vp) == [N[0.1, 0.3], N[0.2, 0.1], N[0.9, 0.2], N[0.4, 0.6]]

    # test support vector of a VPolygon
    p = HPolygon{N}()
    for ci in [c1, c2, c3, c4]
    addconstraint!(p, ci)
    end
    p = tovrep(p)

    # Test Support Vector
    d = N[1., 0.]
    @test σ(d, p) == N[4., 2.]
    d = N[0., 1.]
    @test σ(d, p) == N[2., 4.]
    d = N[-1., 0.]
    @test σ(d, p) == N[-1., 1.]
    d = N[0., -1.]
    @test σ(d, p) == N[0., 0.]

    # test that #83 is fixed
    v = VPolygon([N[2.0, 3.0]])
    @test N[2.0, 3.0] ∈ v
    @test !(N[3.0, 2.0] ∈ v)
    v = VPolygon([N[2.0, 3.0], N[-1.0, -3.4]])
    @test N[-1.0, -3.4] ∈ v

    # an_element function
    v = VPolygon([N[2., 3.]])
    @test an_element(v) ∈ v

    # subset
    p1 = VPolygon([N[0., 0.],N[2., 0.]])
    p2 = VPolygon([N[1., 0.]])
    @test ⊆(p2, p1) && ⊆(p2, p1, true)[1]
    @test ⊆(HPolygon{N}(), p1)
    @test ⊆(p1, BallInf(N[1., 0.], N(1.)))
end
