for N in [Float64, Float32, Rational{Int}]
    # random polygons
    rand(HPolygon)
    rand(HPolygonOpt)
    rand(VPolygon)

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

    # support vector of polygon with no constraints
    @test_throws AssertionError σ(N[0], HPolygon{N}())
    @test_throws AssertionError σ(N[0], HPolygonOpt(HPolygon{N}()))

    # boundedness
    @test isbounded(p)
    @test isbounded(po)

    # isempty
    @test !isempty(p)
    @test !isempty(po)

    # subset
    b1 = Ball1(N[0, 0], N(1))
    b2 = Ball1(N[1, 1], N(4))
    l1 = LinearMap(N[1 0; 0 1], b1)
    l2 = LinearMap(N[1 0; 0 1], b2)
    subset, point = ⊆(p, b1, true)
    @test !subset && point ∈ p && !(point ∈ b1)
    subset, point = ⊆(p, l1, true)
    @test !subset && point ∈ p && !(point ∈ l1)
    subset, point = ⊆(p, b2, true)
    @test subset && ⊆(p, b2) && point == N[]
    subset, point = ⊆(p, l2, true)
    @test subset && ⊆(p, l2) && point == N[]

    # HPolygon/HPolygonOpt tests
    for hp in [p, po]
        # Test Dimension
        @test dim(hp) == 2

        # test support vector, with linear and binary search
        d = N[1, 0]
        @test σ(d, hp) == N[4, 2]
        @test σ(d, hp, linear_search=true) == σ(d, hp, linear_search=false)
        d = N[0, 1]
        @test σ(d, hp) == N[2, 4]
        @test σ(d, hp, linear_search=true) == σ(d, hp, linear_search=false)
        d = N[0, -1]
        @test σ(d, hp) == N[0, 0]
        @test σ(d, hp, linear_search=true) == σ(d, hp, linear_search=false)
        d = N[-1, 0]
        @test σ(d, hp) == N[-1, 1]
        @test σ(d, hp, linear_search=true) == σ(d, hp, linear_search=false)
        d = N[1, -1]
        @test σ(d, hp) == N[4, 2]
        @test σ(d, hp, linear_search=true) == σ(d, hp, linear_search=false)

        # membership
        @test ∈(N[0, 0], hp)
        @test ∈(N[4, 2], hp)
        @test ∈(N[2, 4], hp)
        @test ∈(N[-1, 1], hp)
        @test ∈(N[2, 3], hp)
        @test ∈(N[1, 1], hp)
        @test ∈(N[3, 2], hp)
        @test ∈(N[5 / 4, 7 / 4], hp)
        @test !∈(N[4, 1], hp)
        @test !∈(N[5, 2], hp)
        @test !∈(N[3, 4], hp)
        @test !∈(N[-1, 5], hp)

        # an_element function
        @test an_element(hp) ∈ hp
        hp_shallow = HPolygon{N}()
        @test_throws AssertionError an_element(hp_shallow)
        addconstraint!(hp_shallow, c1)
        @test_throws AssertionError an_element(hp_shallow)

        # hrep conversion
        @test tohrep(hp) == hp

        # constraints list getter
        @test constraints_list(hp) == hp.constraints

        # test constraints list of a VPolygon
        vp = tovrep(hp)
        @test ispermutation(constraints_list(vp), hp.constraints)
    end

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
    @test c == [h1, h4, h2, h5, h3, h6]

    # check that empty polygon (infeasible constraints) has no vertices (#918)
    P = HPolygon([HalfSpace([1.0, 1.0], 0.0), HalfSpace([-1.0, 0.0], -1.0),
        HalfSpace([0.0, -1.0], -1.0)])
    @test vertices_list(P) == Vector{Vector{N}}()

    # is intersection empty
    p3 = convert(HPolygon, Hyperrectangle(low=N[-1, -1], high=N[1, 1]))
    I1 = Interval(N(0.9), N(1.1))
    I2 = Interval(N(0.2), N(0.3))
    I3 = Interval(N(4.0), N(5.0))
    @test !is_intersection_empty(I1 × I2 , p3)
    @test is_intersection_empty(I1 × I3 , p3)

    # Test VRepresentation
    vp = tovrep(p)
    @test N[2, 4] ∈ vertices_list(vp)
    @test N[-1, 1] ∈ vertices_list(vp)
    @test N[0, 0] ∈ vertices_list(vp)
    @test N[4, 2] ∈ vertices_list(vp)
    @test tovrep(vp) == vp

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
    d = N[1, 0]
    @test σ(d, vp) == points[5]
    d = N[0, 1]
    @test σ(d, vp) == points[4]
    d = N[-1, 0]
    @test σ(d, vp) == points[1]
    d = N[0, -1]
    @test σ(d, vp) == points[2]

    # test that #83 is fixed
    v = VPolygon([N[2, 3]])
    @test N[2, 3] ∈ v
    @test !(N[3, 2] ∈ v)
    v = VPolygon([N[2, 3], to_N(N, [-1, -3.4])])
    @test to_N(N, [-1, -3.4]) ∈ v

    # an_element function
    v = VPolygon([N[2, 3]])
    @test an_element(v) ∈ v

    # membership
    point = N[0, 1]
    @test point ∉ VPolygon(Vector{Vector{N}}(), apply_convex_hull=false)
    @test point ∉ VPolygon([N[0, 2]])
    @test point ∈ VPolygon([N[0, 0], N[0, 2]])
    @test point ∉ VPolygon([N[1, 0], N[1, 2]])

    # subset
    p1 = VPolygon([N[0, 0], N[2, 0]])
    p2 = VPolygon([N[1, 0]])
    b = BallInf(N[2, 0], N(1))
    @test ⊆(p2, p1) && ⊆(p2, p1, true)[1]
    @test ⊆(HPolygon{N}(), p1)
    subset, witness = ⊆(p1, b, true)
    @test !⊆(p1, b) && !subset && witness ∈ p1 && witness ∉ b
    subset, witness = ⊆(p2, b, true)
    @test ⊆(p2, b) && subset

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
                @test vi <= v[j]
        end
    end

    # hrep conversion
    v1 = to_N(N, [0.9, 0.2])
    v2 = to_N(N, [0.4, 0.6])
    v3 = to_N(N, [0.2, 0.1])
    v4 = to_N(N, [0.1, 0.3])
    v5 = to_N(N, [0.7, 0.4])
    v6 = to_N(N, [0.8, 0.3])
    v7 = to_N(N, [0.3, 0.55])
    v8 = to_N(N, [0.2, 0.45])
    points5 = [v1, v2, v3, v4, v5, v6, v7, v8]
    for i in [0, 1, 2, 4, 8]
        points = i == 0 ? Vector{Vector{N}}() : points5[1:i]
        vp = VPolygon(points, apply_convex_hull=i > 0)
        h1 = tohrep(vp)
        if i == 0
            @test isempty(h1.constraints)
        elseif i == 1
            @test v1 ∈ h1
        elseif i == 2
            c = h1.constraints[1]
            @test c.a ≈ to_N(N, [0.4, 0.5]) && c.b ≈ to_N(N, (0.46))
            c = h1.constraints[3]
            @test c.a ≈ to_N(N, [-0.4, -0.5]) && c.b ≈ to_N(N, (-0.46))
        elseif i == 4
            @test length(h1.constraints) == 4
            c = h1.constraints[1]
            @test c.a ≈ to_N(N, [0.4, 0.5]) && c.b ≈ to_N(N, (0.46))
            c = h1.constraints[2]
            @test c.a ≈ to_N(N, [-0.3, 0.3]) && c.b ≈ to_N(N, (0.06))
            c = h1.constraints[3]
            @test c.a ≈ to_N(N, [-0.2, -0.1]) && c.b ≈ to_N(N, (-0.05))
            c = h1.constraints[4]
            @test c.a ≈ to_N(N, [0.1, -0.7]) && c.b ≈ to_N(N, (-0.05))
        end

        # check that constraints are sorted correctly
        h2 = HPolygon{N}()
        for c in h1.constraints
            addconstraint!(h2, c)
        end
        @test h1.constraints == h2.constraints
    end
end

function same_constraints(v::Vector{LinearConstraint{N}})::Bool where N<:Real
    c1 = v[1]
    for k = 2:length(v)
        c2 = v[2]
        if c1.a != c2.a || c1.b != c2.b
            return false
        end
    end
    return true
end

for N in [Float64, Float32]
    # test adding constraints, with linear and binary search
    p1 = HPolygon{N}()
    p2 = HPolygon{N}()
    po1 = HPolygonOpt{N}()
    po2 = HPolygonOpt{N}()
    n = 10
    for i in 1:n
        constraint = LinearConstraint(rand(N, 2), rand(N))
        addconstraint!(p1, constraint, linear_search=true)
        addconstraint!(p2, constraint, linear_search=i<=2)
        addconstraint!(po1, constraint, linear_search=true)
        addconstraint!(po2, constraint, linear_search=i<=2)
    end
    for i in 1:n
        @test same_constraints([p1.constraints[i], p2.constraints[i],
                                po1.constraints[i], po2.constraints[i]])
    end
end

# default Float64 constructors
@test HPolygon() isa HPolygon{Float64}
@test HPolygonOpt() isa HPolygonOpt{Float64}
@test VPolygon() isa VPolygon{Float64}
