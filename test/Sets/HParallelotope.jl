for N in [Float32, Float64, Rational{Int}]
    # See [DreossiDP17; Example 6](@citet).
    D = N[-1 0 0; -1 -1 0; 0 0 -1]
    c = N[-0.80, -0.95, 0, 0.85, 1, 0]
    P = HParallelotope(D, c)

    # illegal input: empty set
    D2 = N[1 1; 0 1]
    c2 = N[5 // 2, 2, -4, -1]
    @test_throws ArgumentError HParallelotope(D2, c2)
    P2 = HParallelotope(D2, c2; check_consistency=false)
    @test isempty(convert(HPolyhedron, P2))

    # illegal input: unbounded set
    D2 = N[1 1; 1 1]
    c2 = N[1, 1, -1, -1]
    @test_throws ArgumentError HParallelotope(D2, c2)
    P2 = HParallelotope(D2, c2; check_consistency=false)
    @test !isbounded(convert(HPolyhedron, P2))

    # test getter functions
    @test directions(P) == D
    @test offset(P) == c
    @test dim(P) == 3

    # ispolyhedral
    @test ispolyhedral(P)

    # check constructor's size check
    @test_throws AssertionError HParallelotope(D, vcat(c, c))

    # test computation of the list of constraints
    @test constraints_list(P) == [HalfSpace(D[1, :], c[1]),
                                  HalfSpace(D[2, :], c[2]),
                                  HalfSpace(D[3, :], c[3]),
                                  HalfSpace(-D[1, :], c[4]),
                                  HalfSpace(-D[2, :], c[5]),
                                  HalfSpace(-D[3, :], c[6])]
    P3 = HParallelotope(zeros(N, 0, 0), zeros(N, 0); check_consistency=false)
    @test isempty(constraints_list(P3))

    P = HParallelotope(N[1 0; 0 1], N[1, 1, 1, 1])

    # test center
    @test center(P) == N[0, 0]

    # test vertices functions
    @test base_vertex(P) == N[-1, -1]
    @test extremal_vertices(P) == [N[1, -1], N[-1, 1]]

    # test generators getters
    @test genmat(P) == N[1 0; 0 1]

    # random parallelotope
    rand(HParallelotope)

    # emptiness
    @test !isempty(P)

    # volume
    @test volume(P) == N(4)

    # vertices list
    @test ispermutation(vertices_list(P),
                        [N[1, 1], N[1, -1], N[-1, 1], N[-1, -1]])
end

for N in [Float32, Float64]
    # conversion from zonotope
    Z = Zonotope(N[0, 0], N[1 0; 0 1])
    @test convert(HParallelotope, Z) ==
          HParallelotope(N[0 -1; 1 0], N[1, 1, 1, 1])
    Z = Zonotope(N[0, 0], N[1 0 1; 0 1 1])
    @test_throws AssertionError convert(HParallelotope, Z)
end

for N in [Float64]
    P = HParallelotope(N[1 0; 0 1], N[1, 1, 1, 1])
    @test collect(generators(P)) == [N[1, 0], N[0, 1]]
end

@test !isoperationtype(HParallelotope)
