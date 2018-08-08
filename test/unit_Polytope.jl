for N in [Float64, Rational{Int}, Float32]
    # -----
    # H-rep
    # -----

    # constructor from matrix and vector
    A = [N(1.) N(2.); N(-1.) N(1.)]
    b = [N(1.), N(2.)]
    p = HPolytope(A, b)
    c = p.constraints
    @test c isa Vector{LinearConstraint{N}}
    @test c[1].a == N[1.0, 2.0] && c[1].b == N(1.0)
    @test c[2].a == N[-1.0, 1.0] && c[2].b == N(2.0)

    # convert back to matrix and vector
    A2, b2 = tosimplehrep(p)
    @test A == A2 && b == b2

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
    @test ∈(N[5. / 4., 7. / 4.], p)
    @test !∈(N[4., 1.], p)

    # singleton list (only available with Polyhedra library)
    if !isdefined(@__MODULE__, :Polyhedra)
        @test_throws MethodError singleton_list(p)
    end

    # -----
    # V-rep
    # -----

    # constructor from a VPolygon
    polygon = VPolygon([N[0., 0.], N[1., 0.], N[0., 1.]])
    p = VPolytope(polygon)
    @test vertices_list(polygon) == vertices_list(p)

    # dim
    @test dim(p) == 2

    # support vector
    d = N[1., 0.]
    @test_throws ErrorException σ(d, p, algorithm="xyz")
    if !isdefined(@__MODULE__, :Polyhedra)
        @test_throws AssertionError σ(d, p)
    end

    # vertices_list function
    @test vertices_list(p) == p.vertices
end

# default Float64 constructors
@test HPolytope() isa HPolytope{Float64}
@test VPolytope() isa VPolytope{Float64}
