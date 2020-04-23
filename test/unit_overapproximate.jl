using LazySets.Approximations: project,
                               get_linear_coeffs

for N in [Float64, Rational{Int}, Float32]
    c = N[0, 0]
    b = Ball1(c, N(1))

    # overapproximating a set of type T1 with an unsupported type T2 is the
    # identity if T1 = T2
    @test_throws MethodError overapproximate(ZeroSet{N}(2), EmptySet)
    e = EmptySet{N}(2)
    @test overapproximate(e, EmptySet) == e
    # the third argument is ignored
    oa = overapproximate(b, Ball1, Hyperrectangle)
    @test oa isa Ball1

    # HPolygon approximation with box directions
    p = overapproximate(b, HPolygon)
    for d in [N[1, 0], N[-1, 0]]
        @test σ(d, p)[1] ≈ σ(d, b)[1]
    end
    for d in [N[0, 1], N[0, -1]]
        @test σ(d, p)[2] ≈ σ(d, b)[2]
    end

    # Hyperrectangle approximation
    p = overapproximate(b, Hyperrectangle)
    for d in [N[1, 0], N[-1, 0]]
        @test σ(d, p)[1] ≈ σ(d, b)[1]
    end
    for d in [N[0, 1], N[0, -1]]
        @test σ(d, p)[2] ≈ σ(d, b)[2]
    end
    @test p.center ≈ c
    @test p.radius ≈ N[1, 1]

    # Interval approximation
    b = Ball1(N[0], N(1))
    p1 = overapproximate(b, LazySets.Interval)
    cap = Intersection(b, BallInf(N[1], N(1)))
    p2 = overapproximate(cap, LazySets.Interval)
    caparray = IntersectionArray([b, BallInf(N[1], N(1)), BallInf(N[1//2], N(1//4))])
    p3 = overapproximate(caparray, LazySets.Interval)
    for d in [N[1], N[-1]]
        for (o, p) in [(b, cap, caparray), (p1, p2, p3)]
            @test σ(d, o)[1] ≈ d[1]
        end
    end

    # approximation with an axis-aligned hyperrectangle
    Z = Zonotope(N[-1.0, -1.0], N[-1/2 0; -1.2 -1])
    Zoa = overapproximate(Z, Hyperrectangle) # faster o.a.
    Zba = box_approximation(Z) # default o.a. implementation that uses supp function
    @test Zoa.center ≈ Zba.center && Zoa.radius ≈ Zba.radius

    # same but using cartesian product
    h1 = Hyperrectangle(N[1/2],  N[1/2])
    h2 = Hyperrectangle(N[2.5, 4.5],  N[1/2, 1/2])
    H = overapproximate(h1 × h2, Hyperrectangle) # defaults to convert method
    @test low(H) == N[0, 2, 4] && high(H) == N[1, 3, 5]

    # overapproximation of the lazy linear map of a hyperrectangular set
    H = Hyperrectangle(N[0, 0], N[1/2, 1])
    M = Diagonal(N[2, 2])
    OA = overapproximate(M*H, Hyperrectangle)
    @test OA isa Hyperrectangle && OA.center == N[0, 0] && OA.radius == N[1, 2]

    #overapproximation of Minkowski sum of linear maps for each block in the row block
    i1 = Interval(N[0, 1])
    h = Hyperrectangle(low=N[3, 4], high=N[5, 7])
    M = N[1 2 3; 4 5 6; 7 8 9]
    cpa = CartesianProductArray([i1, h])
    lm = M * cpa

    oa = overapproximate(lm, Hyperrectangle)
    oa_box = overapproximate(lm, Approximations.BoxDirections)
    d_oa_d_hp = overapproximate(lm, CartesianProductArray{N, Hyperrectangle{N}})
    d_oa_d_box = overapproximate(lm, CartesianProductArray, Approximations.BoxDirections)
    oa_d_hp = overapproximate(d_oa_d_hp)

    @test oa == oa_d_hp

    for (oax, set_type) in [(d_oa_d_hp, Hyperrectangle), (d_oa_d_box, HPolytope)]
        @test oax isa CartesianProductArray
        arr = oax.array
        @test length(arr) == 2 && dim(arr[1]) == 1 && dim(arr[2]) == 2
        @test all(X -> X isa set_type, arr)
    end

    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    i3 = Interval(N[1, 4])
    cpa = CartesianProductArray([i1, i2, i3])
    M = N[1 2 0; 0 1 0; 0 1 1]
    lm = M * cpa
    d_oa = overapproximate(lm, CartesianProductArray{N, Interval{N}})
    oa = overapproximate(lm)
    @test overapproximate(d_oa) == oa
    @test typeof(d_oa) == CartesianProductArray{N, Interval{N}}

    # interval linear map of zonotope
    M = IntervalMatrix([(N(0) ± N(1)) (N(0) ± N(0)); (N(0) ± N(0)) (N(0) ± N(1))])
    Z = Zonotope(N[1, 1], N[1 1; 1 -1])
    lm = LinearMap(M, Z)
    Zo = overapproximate(lm, Zonotope)
    @test box_approximation(Zo) == Hyperrectangle(N[0, 0], N[3, 3])

    # rectification
    r = Rectification(EmptySet{N}(2))
    @test overapproximate(r, Hyperrectangle) isa EmptySet{N}
    r = Rectification(Ball1(N[0, 0], N(1)))
    @test overapproximate(r, Hyperrectangle) ==
        Hyperrectangle(low=N[0, 0], high=N[1, 1])

    # HPolytope error messages
    @test_throws ArgumentError overapproximate(Singleton(N[0, 0]), HPolytope)
    @test_throws ArgumentError overapproximate(Singleton(N[0]), HPolytope)

    # Hyperrectangle of AffineMap
    A = N[1 2; 1 3; 1 4]; X = Hyperrectangle(N[0, 1], N[0, 1]); b = N[1, 2, 3];
    am = AffineMap(A, X, b)
    @test overapproximate(am, Hyperrectangle) == Hyperrectangle(N[3, 5, 7], N[2, 3, 4])
end

# tests that do not work with Rational{Int}
for N in [Float64, Float32]
    # useful for benchmarking overapproximate and LinearMap's support vector
    # (see #290)
    function overapproximate_lmap(n)
        B = BallInf(ones(N, n), N(2))
        π = sparse(N[1, 2], N[1, 2], ones(N, 2), 2, n)
        return Approximations.overapproximate(π*B)
    end
    o = overapproximate_lmap(50)
    @test o.center == N[1, 1] && o.radius == N[2, 2]

    # Approximation of a 2D centered unit ball in norm 1
    # All vertices v should be like this:
    # ‖v‖ >= 1 and ‖v‖ <= 1+ε
    # Where ε is the given error bound
    b = Ball1(N[0, 0], N(1))
    ε = N(0.01)
    p = tovrep(overapproximate(b, ε))
    for v in vertices_list(p)
    @test norm(v) >= N(1)
    @test norm(v) <= N(1 + ε)
    end

    # Check that there are no redundant constraints for a ballinf
    b = BallInf(N[0.5, 0.5], N(0.1))
    lcl = overapproximate(b, N(0.001)).constraints
    @test length(lcl) == 4
    @test lcl[1].a == N[1, 0]
    @test lcl[1].b == N(0.6)
    @test lcl[2].a == N[0, 1]
    @test lcl[2].b == N(0.6)
    @test lcl[3].a == N[-1, 0]
    @test lcl[3].b == N(-0.4)
    @test lcl[4].a == N[0, -1]
    @test lcl[4].b == N(-0.4)

    # Check that there are no redundant constraints for a HPolygon (octagon)
    p = HPolygon{N}()
    addconstraint!(p, LinearConstraint(N[1, 0], N(1)))
    addconstraint!(p, LinearConstraint(N[0, 1], N(1)))
    addconstraint!(p, LinearConstraint(N[-1, 0], N(1)))
    addconstraint!(p, LinearConstraint(N[0, -1], N(1)))
    addconstraint!(p, LinearConstraint(N[sqrt(2)/2, sqrt(2)/2], N(1)))
    addconstraint!(p, LinearConstraint(N[-sqrt(2)/2, sqrt(2)/2], N(1)))
    addconstraint!(p, LinearConstraint(N[sqrt(2)/2, -sqrt(2)/2], N(1)))
    addconstraint!(p, LinearConstraint(N[-sqrt(2)/2, -sqrt(2)/2], N(1)))
    lcl = overapproximate(p, N(.001)).constraints
    @test length(lcl) == 8
    @test lcl[1].a ≈ N[1, 0]
    @test lcl[1].b ≈ N(1)
    @test lcl[2].a ≈ N[sqrt(2)/2, sqrt(2)/2]
    @test lcl[2].b ≈ N(1)
    @test lcl[3].a ≈ N[0, 1]
    @test lcl[3].b ≈ N(1)
    @test lcl[4].a ≈ N[-sqrt(2)/2, sqrt(2)/2]
    @test lcl[4].b ≈ N(1)
    @test lcl[5].a ≈ N[-1, 0]
    @test lcl[5].b ≈ N(1)
    @test lcl[6].a ≈ N[-sqrt(2)/2, -sqrt(2)/2]
    @test lcl[6].b ≈ N(1)
    @test lcl[7].a ≈ N[0, -1]
    @test lcl[7].b ≈ N(1)
    @test lcl[8].a ≈ N[sqrt(2)/2, -sqrt(2)/2]
    @test lcl[8].b ≈ N(1)

    # Zonotope approximation of convex hull of two zonotopes

    # same order (mean method only)
    Z1 = Zonotope(zeros(N, 2), hcat(N[1, 0]))
    Z2 = Zonotope(zeros(N, 2), hcat(N[0, 1]))
    Zch = overapproximate(ConvexHull(Z1, Z2), Zonotope; algorithm="mean")
    # the result is a diamond (tight)
    @test Zch == Zonotope(N[0, 0], N[0.5 0.5; 0.5 -0.5])

    # different order (mean method only)
    Z1 = Zonotope(zeros(N, 2), hcat(N[1, 0]))
    Z2 = Zonotope(N[0, 0], N[1 0; 0 1])
    # the result is a box (tight)
    Zch = overapproximate(ConvexHull(Z1, Z2), Zonotope; algorithm="mean")
    @test Zch == Zonotope(N[0, 0], N[1 0; 0 1])

    # comparison of different methods
    Z1 = Zonotope(ones(N, 2), [N[1, 0], N[0, 1], N[1, 1]])
    Z2 = Zonotope(-ones(N, 2), [N[0.5, 1], N[-0.1, 0.9], N[1, 4]])
    Y = ConvexHull(Z1, Z2)
    # overapproximate with a polygon
    Y_polygon = overapproximate(Y, N(1e-3))
    # overapproximate with a zonotope (different methods)
    overapproximate(Y, Zonotope) # default
    Y_zonotope_1 = overapproximate(Y, Zonotope; algorithm="mean")
    Y_zonotope_2 = overapproximate(Y, Zonotope; algorithm="join")
    @test Y_polygon ⊆ Y_zonotope_1
    @test Y_polygon ⊆ Y_zonotope_2
    # check property of 'join' algorithm: box approximation equality
    X1 = box_approximation(Y)
    X2 = box_approximation(Y_zonotope_2)
    @test center(X1) ≈ center(X2) &&
          radius_hyperrectangle(X1) ≈ radius_hyperrectangle(X2)

    # example from paper
    Z1 = Zonotope(N[3, 0], N[1 2; 1 1])
    Z2 = Zonotope(N[1, 0], N[-2 1; 1 1])
    @test overapproximate(ConvexHull(Z1, Z2), Zonotope; algorithm="join") ==
          Zonotope(N[2, 0], N[0 1 3; 1 1 0])

    # different orders
    Z1 = Zonotope(N[3, 0], N[1 2 1; 1 1 2])
    Z2 = Zonotope(N[1, 0], N[-2 1; 1 1])
    Y = ConvexHull(Z1, Z2)
    Y_polygon = overapproximate(Y, N(1e-3))
    Y_zonotope_1 = overapproximate(Y, Zonotope; algorithm="mean")
    Y_zonotope_2 = overapproximate(Y, Zonotope; algorithm="join")
    @test Y_polygon ⊆ Y_zonotope_1
    @test Y_polygon ⊆ Y_zonotope_2

    # ResetMap and CPA
    X = CartesianProductArray([Ball1(N[0, 0], N(1)), Ball1(N[0, 0], N(2)),
                               Ball1(N[0, 0], N(3))])
    resets = Dict(3 => N(0), 6 => N(0))
    rm = ResetMap(X, resets)
    Y = overapproximate(rm, CartesianProductArray, Hyperrectangle)
    @test array(Y) == [Hyperrectangle(N[0, 0], N[1, 1]),
        Hyperrectangle(N[0, 0], N[0, 2]), Hyperrectangle(N[0, 0], N[3, 0])]
end

for N in [Float64]
    i1 = Interval(N[0, 1])
    h = Hyperrectangle(low=N[3, 4], high=N[5, 7])
    M = N[1 2 3; 4 5 6; 7 8 9]
    cpa = CartesianProductArray([i1, h])
    lm = M * cpa
    d_oa_d_box = overapproximate(lm, CartesianProductArray, Approximations.BoxDirections)
    oa_d_box = overapproximate(d_oa_d_box, Approximations.BoxDirections)
    oa_box = overapproximate(lm, Approximations.BoxDirections)
    @test oa_box == oa_d_box

    # intersection of two polyhedra
    P = HPolyhedron([HalfSpace(N[1, 0], N(1)), HalfSpace(N[0, 1], N(1))])
    Q = HPolyhedron([HalfSpace(N[-1, 0], N(1)), HalfSpace(N[0, -1], N(1))])
    oa = overapproximate(P ∩ Q, BoxDirections)
    B = BallInf(N[0, 0], N(1))
    @test B ⊆ oa && oa ⊆ B

    # decomposed linear map approximation
    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    M = N[1 2; 0 1]
    cpa = CartesianProductArray([i1, i2])
    lm = M * cpa
    d_oa = overapproximate(lm, CartesianProductArray{N, Interval{N}})
    oa = overapproximate(lm, OctDirections)
    @test oa ⊆ d_oa

    # decomposed intersection between Cartesian product array and polyhedron
    h1 = Hyperrectangle(low=N[3, 4], high=N[5, 7])
    h2 = Hyperrectangle(low=N[5, 5], high=N[6, 8])
    cpa1 = CartesianProductArray([i1, i2, h1])
    o_cpa1 = overapproximate(cpa1)
    cpa2 = CartesianProductArray([i1, i2, h2])
    o_cpa2 = overapproximate(cpa2)
    G = HalfSpace(N[1, 0, 0, 0], N(1))
    G_3 = HalfSpace(N[0, 0, 1, 0], N(1))
    G_3_neg = HalfSpace(N[0, 0, -1, 0], N(0))
    G_comb = HalfSpace(N[1, 1, 0, 0], N(2.5))

    d_int_g = Intersection(cpa1, G)
    d_int_g_3 = Intersection(cpa1, G_3)
    d_int_g_3_neg = Intersection(cpa1, G_3_neg)
    d_int_cpa = intersection(cpa1, cpa2)
    l_int_cpa = intersection(o_cpa1, o_cpa2)
    o_d_int_g = overapproximate(d_int_g, CartesianProductArray, Hyperrectangle)
    o_d_int_g_3 = overapproximate(d_int_g_3, CartesianProductArray, Hyperrectangle)
    o_d_int_g_3_neg = overapproximate(d_int_g_3_neg, CartesianProductArray, Hyperrectangle)

    @test overapproximate(o_d_int_g) ≈ overapproximate(d_int_g)
    @test isempty(d_int_g_3) && o_d_int_g_3 == EmptySet{N}(4)
    @test overapproximate(o_d_int_g_3_neg) ≈ overapproximate(d_int_g_3_neg)
    @test overapproximate(d_int_cpa, Hyperrectangle) == l_int_cpa
    @test all([X isa CartesianProductArray for X in [d_int_cpa, o_d_int_g, o_d_int_g_3_neg]])

    int_g_comb = Intersection(cpa1, G_comb)
    o_d_int_g_comb = overapproximate(int_g_comb, CartesianProductArray, Hyperrectangle)
    o_bd_int_g_comb = overapproximate(int_g_comb, BoxDirections)
    @test overapproximate(o_d_int_g_comb) ≈ overapproximate(o_bd_int_g_comb)
    projection = project(o_bd_int_g_comb, [1, 2], BoxDirections)
    @test vertices_list(projection) ≈ [[0.5, 2.5], [0., 2.5], [0., 2], [0.5, 2.]]
    projection = project(overapproximate(int_g_comb, OctDirections), [1, 2], OctDirections)
    @test vertices_list(projection) ≈ [[0., 2.5], [0., 2.], [0.5, 2.]]

    # Zonotope approximation
    Z1 = Zonotope(ones(N, 2), [N[1, 0], N[0, 1], N[1, 1]])
    Z2 = Zonotope(-ones(N, 2), [N[0.5, 1], N[-0.1, 0.9], N[1, 4]])
    Y = ConvexHull(Z1, Z2)
    Y_polygon = overapproximate(Y, N(1e-3)) # overapproximate with a polygon
    Y_zonotope = overapproximate(Y, Zonotope) # overapproximate with a zonotope
    @test Y_polygon ⊆ Y_zonotope

    # Zonotope approximation with fixed generator directions
    X = Ball1(zeros(N, 2), N(1))
    Y = overapproximate(X, Zonotope, BoxDirections{N}(2))
    @test Y == Zonotope(N[0, 0], N[1 0; 0 1])
    Y = overapproximate(X, Zonotope, OctDirections{N}(2))
    @test remove_zero_generators(Y) == Zonotope(N[0, 0], N[1//2 1//2; 1//2 -1//2])
    X = Ball1(zeros(N, 3), N(1))
    Y = overapproximate(X, Zonotope, BoxDirections{N}(3))
    @test Y == Zonotope(N[0, 0, 0], N[1 0 0; 0 1 0; 0 0 1])

    # decomposed linear map approximation
    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    M = N[1 2; 0 1]
    cpa = CartesianProductArray([i1, i2])
    lm = M * cpa
    d_oa = overapproximate(lm, CartesianProductArray{N, Interval{N}})
    oa = overapproximate(lm, OctDirections)
    @test oa ⊆ d_oa

    # =======================================
    # Zonotope overapprox. of a Taylor model
    # =======================================
    x₁, x₂, x₃ = set_variables(N, ["x₁", "x₂", "x₃"], order=5)
    Dx₁ = IA.Interval(N(1.0), N(3.0))
    Dx₂ = IA.Interval(N(-1.0), N(1.0))
    Dx₃ = IA.Interval(N(-1.0), N(0.0))
    D = Dx₁ × Dx₂ × Dx₃   # domain
    x0 = IntervalBox(IA.mid.(D)...)
    I = IA.Interval(N(0.0), N(0.0)) # interval remainder
    p₁ = 1 + x₁ - x₂
    p₂ = x₃ - x₁
    vTM = [TaylorModels.TaylorModelN(pi, I, x0, D) for pi in [p₁, p₂]]
    Z1 = overapproximate(vTM, Zonotope)
    @test center(Z1) == N[3, -2.5]
    @test Matrix(genmat(Z1)) == N[1 -1 0; -1 0 0.5]

    # auxiliary function to get the linear coefficients
    t = TaylorModels.Taylor1(0) # t.order is 0
    @test get_linear_coeffs(t) == N[0]
    p = x₁ + 2x₂ - 3x₃
    @test get_linear_coeffs(p) == N[1, 2, -3]
    y = set_variables("y", numvars=2, order=1)
    p = zero(y[1])
    @test get_linear_coeffs(p) == N[0, 0]
end
