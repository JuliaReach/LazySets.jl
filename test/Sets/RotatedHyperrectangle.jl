for N in [Float64, Rational{Int}, Float32]
    c = N[1, 1]
    r = N[1, 1]
    H = Hyperrectangle(c, r)
    M = N[2 1; 1 -1]
    R = RotatedHyperrectangle(M, H)

    # dimension
    @test dim(R) == 2
    # center
    @test center(R) == N[3, 0]

    # generators
    gens = collect(generators(R))
    @test ngens(R) == 2
    @test length(gens) == 2 &&
        (N[2, 1] ∈ gens || N[-2, -1] ∈ gens) &&
        (N[1, -1] ∈ gens || N[-1, 1] ∈ gens)
    @test genmat(R) ∈ [N[2 1; 1 -1], N[-2 1; -1 -1], N[2 -1; 1 1], N[-2 -1; -1 1]]

    # support vector and support function
    d = N[1, 0]
    @test σ(d, R) == N[6, 0]
    @test ρ(d, R) == N(6)
    d = N[-1, 0]
    @test σ(d, R) == N[0, 0]
    @test ρ(d, R) == N(0)
    d = N[0, 1]
    @test σ(d, R) == N[4, 2]
    @test ρ(d, R) == N(2)
    d = N[0, -1]
    @test σ(d, R) == N[2, -2]
    @test ρ(d, R) == N(2)

    # boundedness
    @test isbounded(R)

    # is_polyhedral
    @test is_polyhedral(R)

    # isempty
    @test !isempty(R)

    # universality
    answer, w = isuniversal(R, true)
    @test !isuniversal(R) && !answer && w ∉ R

    # membership
    @test N[4, 1] ∈ R
    @test N[6, 1] ∉ R

    # an_element function
    @test an_element(R) ∈ R

    # linear map
    R2 = linear_map(M, R)
    @test R2 == RotatedHyperrectangle(M^2, H)
    R2 = linear_map(N[1 1;], R)
    @test R2 == RotatedHyperrectangle(N[3 0;], H)

    # list of vertices
    @test ispermutation(vertices_list(R),
                        [N[6, 0], N[2, -2], N[0, 0], N[4, 2]])

    if test_suite_polyhedra
        # list of constraints (only the length, see #2626)
        @test length(constraints_list(R)) == 4
    end
end
