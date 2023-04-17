for N in [Float64, Float32, Rational{Int}]
    # constructor
    vertices = [N[1, 0, -1 / sqrt(2)], N[-1, 0, -1 / sqrt(2)], N[0, 1, 1 / sqrt(2)],
                N[0, -1, 1 / sqrt(2)]]
    T = Tetrahedron(vertices)
    @test eltype(T) == N
    @test dim(T) == 3

    vmat = reduce(hcat, vertices)
    Tmat = Tetrahedron(vmat)
    @test Tmat == T

    # suport vector
    @test σ(ones(3), T) == N[0, 1, 1 / sqrt(2)]

    # membership
    @test zeros(N, 3) ∈ T

    # is_polyhedral
    @test is_polyhedral(T)

    # LazySets#3303
    T = Tetrahedron([N[0, 0, 0], N[0, 1, 0], N[1, 0, 0], N[1, 0, 1]])
    @test zeros(N, 3) ∈ T
end
