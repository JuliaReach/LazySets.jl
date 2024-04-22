# check solver handling
for N in [Float64, Rational{Int}, Int]
    LazySets.default_lp_solver(N)
end

# check that the default model and explicit solver specification work for linprog
for N in [Float64]
    p = HPolyhedron{N}()
    c1 = LinearConstraint(N[2, 2], N(12))
    c2 = LinearConstraint(N[-3, 3], N(6))
    c3 = LinearConstraint(N[-1, -1], N(0))
    c4 = LinearConstraint(N[2, -4], N(0))
    addconstraint!(p, c3)
    addconstraint!(p, c1)
    addconstraint!(p, c4)
    addconstraint!(p, c2)

    P = HPolyhedron([HalfSpace(N[3 // 50, -4 // 10], N(1)),
                     HalfSpace(N[-1 // 50, 1 // 10], N(-1))])

    ####################
    # with default model
    ####################

    # support vector
    d = N[1, 0]
    @test σ(d, p) == N[4, 2]
    d = N[0, 1]
    @test σ(d, p) == N[2, 4]
    d = N[-1, 0]
    @test σ(d, p) == N[-1, 1]
    d = N[0, -1]
    @test σ(d, p) == N[0, 0]

    # an_element
    @test an_element(P) ∈ P

    ######################
    # with explicit solver
    ######################
    solver = optimizer_with_attributes(() -> GLPK.Optimizer(; method=GLPK.SIMPLEX))

    # support vector
    d = N[1, 0]
    @test σ(d, p; solver=solver) == N[4, 2]
    d = N[0, 1]
    @test σ(d, p; solver=solver) == N[2, 4]
    d = N[-1, 0]
    @test σ(d, p; solver=solver) == N[-1, 1]
    d = N[0, -1]
    @test σ(d, p; solver=solver) == N[0, 0]

    # an_element
    @test an_element(P; solver=solver) ∈ P
end
