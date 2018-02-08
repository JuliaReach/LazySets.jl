for N in [Float64, Float32, Rational{Int}]
    # Cartesian Product of a centered 1D BallInf and a centered 2D BallInf
    # Here a 3D BallInf
    b1 = BallInf(N[0.], N(1.))
    b2 = BallInf(N[0., 0.], N(1.))
    # Test Construction
    c1 = CartesianProduct(b1, b2)
    @test c1.X == b1
    @test c1.Y == b2
    # Test Dimension
    @test dim(c1) == 3
    # Test Support Vector
    d = N[1., 1., 1.]
    @test σ(d, c1) == N[1., 1., 1.]
    d = N[1., 1., -1.]
    @test σ(d, c1) == N[1., 1., -1.]
    d = N[1., -1., 1.]
    @test σ(d, c1) == N[1., -1., 1.]
    d = N[1., -1., -1.]
    @test σ(d, c1) == N[1., -1., -1.]
    d = N[-1., 1., 1.]
    @test σ(d, c1) == N[-1., 1., 1.]
    d = N[-1., 1., -1.]
    @test σ(d, c1) == N[-1., 1., -1.]
    d = N[-1., -1., 1.]
    @test σ(d, c1) == N[-1., -1., 1.]
    d = N[-1., -1., -1.]
    @test σ(d, c1) == N[-1., -1., -1.]

    # Cartesian Product of a not-centered 1D BallInf and a not-centered 2D BallInf
    # Here a Hyperrectangle where c = [1., -3., 4.] and r = [3., 2., 2.]
    b1 = BallInf(N[1.], N(3.))
    b2 = BallInf(N[-3., 4.], N(2.))
    # Test Construction
    c = CartesianProduct(b1, b2)
    # Test Dimension
    @test dim(c) == 3
    # Test Support Vector
    d = N[1., 1., 1.]
    @test σ(d, c) == N[4., -1., 6.]
    d = N[1., 1., -1.]
    @test σ(d, c) == N[4., -1., 2.]
    d = N[1., -1., 1.]
    @test σ(d, c) == N[4., -5., 6.]
    d = N[1., -1., -1.]
    @test σ(d, c) == N[4., -5., 2.]
    d = N[-1., 1., 1.]
    @test σ(d, c) == N[-2., -1., 6.]
    d = N[-1., 1., -1.]
    @test σ(d, c) == N[-2., -1., 2.]
    d = N[-1., -1., 1.]
    @test σ(d, c) == N[-2., -5., 6.]
    d = N[-1., -1., -1.]
    @test σ(d, c) == N[-2., -5., 2.]

    # Test Cartesian Product with EmptySet
    s = Singleton(N[1.])
    E = EmptySet{N}()
    cs1 = E * s
    cs2 = s * E
    @test cs1 isa EmptySet
    @test cs2 isa EmptySet

    # Test Cartesian Product of an array
    # 0-elements
    as = LazySet{N}[]
    cs = CartesianProduct(as)
    @test cs isa EmptySet
    # 1-element
    as = [Singleton(N[1.])]
    cs = CartesianProduct(as)
    @test cs.element == N[1.]
    # 3-elements
    as = [Singleton(N[1.]), Singleton(N[2.]), Singleton(N[3.])]
    cs = CartesianProduct(as)
    @test cs.X.element == N[1.]
    @test cs.Y.X.element == N[2.]
    @test cs.Y.Y.element == N[3.]

    # Test containment with respect to CartesianProduct
    p1 = HPolygon{N}()
    addconstraint!(p1, LinearConstraint(N[2., 2.], N(12.)))
    addconstraint!(p1, LinearConstraint(N[-3., 3.], N(6.)))
    addconstraint!(p1, LinearConstraint(N[-1., -1.], N(0.)))
    addconstraint!(p1, LinearConstraint(N[2., -4.], N(0.)))
    p2 = HPolygon{N}()
    addconstraint!(p2, LinearConstraint(N[1., 0.], N(1.)))
    addconstraint!(p2, LinearConstraint(N[-1., 0.], N(1.)))
    addconstraint!(p2, LinearConstraint(N[0., 1.], N(1.)))
    addconstraint!(p2, LinearConstraint(N[0., -1.], N(1.)))
    cp = CartesianProduct(p1, p2)

    @test ∈(N[0., 0., 0., 0.], cp)
    @test ∈(N[4., 2., 1., 0.], cp)
    @test ∈(N[2., 4., -1., 0.], cp)
    @test ∈(N[-1., 1., .5, .7], cp)
    @test ∈(N[2., 3., -.8, .9], cp)
    @test ∈(N[1., 1., -1., 0.], cp)
    @test ∈(N[3., 2., 0., 1.], cp)
    @test ∈(N[5./4., 7./4., 1., 1.], cp)
    @test !∈(N[4., 1., 0., 0.], cp)
    @test !∈(N[5., 2., 0., 0.], cp)
    @test !∈(N[3., 4., 0., 0.], cp)
    @test !∈(N[-1., 5., 0., 0.], cp)
    @test !∈(N[4., 2., 3., 1.], cp)
    @test !∈(N[2., 3., 3., 1.], cp)
    @test !∈(N[0., 0., 3., 1.], cp)
    @test !∈(N[1., 1., 3., 1.], cp)

    # array getter
    v = Vector{N}(0)
    @test array(CartesianProductArray()) == v
end
