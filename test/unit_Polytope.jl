for N in [Float64, Rational{Int}, Float32]
    # 2D polytope
    p = HPolytope{N}()
    c1 = LinearConstraint(N[2., 2.], N(12.))
    c2 = LinearConstraint(N[-3., 3.], N(6.))
    c3 = LinearConstraint(N[-1., -1.], N(0.))
    c4 = LinearConstraint(N[2., -4.], N(0.))
    addconstraint!(p, c3)
    addconstraint!(p, c1)
    addconstraint!(p, c4)
    addconstraint!(p, c2)

    # support vector
    d = N[1., 0.]
    @test σ(d, p) == N[4., 2.]
    d = N[0., 1.]
    @test σ(d, p) == N[2., 4.]
    d = N[-1., 0.]
    @test σ(d, p) == N[-1., 1.]
    d = N[0., -1.]
    @test σ(d, p) == N[0., 0.]

    # membership
    @test ∈(N[5./4., 7./4.], p)
    @test !∈(N[4., 1.], p)

    # singleton list (only available with Polyhedra library)
    @test_throws MethodError singleton_list(p)
end
