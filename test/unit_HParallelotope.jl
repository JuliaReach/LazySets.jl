for N in [Float32, Float64, Rational{Int}]
    # Example 6 from [1].
    #
    # [1] Tommaso Dreossi, Thao Dang, and Carla Piazza. *Reachability computation for polynomial dynamical systems.*
    #      Formal Methods in System Design 50.1 (2017): 1-38.
    D = N[-1 0 0; -1 -1 0; 0 0 -1]
    c = N[-0.80, -0.95, 0.0, 0.85, 1.0, 0.0]
    P = HParallelotope(D, c)

    # test getter functions
    @test directions(P) == D
    @test offset(P) == c
    @test dim(P) == 3

    # check constructor's size check
    @test_throws AssertionError HParallelotope(D, vcat(c, c))

    # test computation of the list of constraints
    @test constraints_list(P) == [HalfSpace(D[1, :], c[1]),
                                  HalfSpace(D[2, :], c[2]),
                                  HalfSpace(D[3, :], c[3]),
                                  HalfSpace(-D[1, :], c[4]),
                                  HalfSpace(-D[2, :], c[5]),
                                  HalfSpace(-D[3, :], c[6])]

    P = HParallelotope(N[1 0; 0 1], N[1, 1, 1, 1])

    # test center
    @test center(P) == N[0, 0]

    # test vertices functions
    @test base_vertex(P) == N[-1, -1]
    @test extremal_vertices(P) == [N[1, -1], N[-1, 1]]

    # test generators getters
    @test genmat(P) == N[1 0; 0 1]
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
