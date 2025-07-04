using LazySets, Test

for N in [Float64, Float32, Rational{Int}]
    # constructor
    vertices = [N[1, 0, -1 / sqrt(2)], N[-1, 0, -1 / sqrt(2)], N[0, 1, 1 / sqrt(2)],
                N[0, -1, 1 / sqrt(2)]]
    T = Tetrahedron(vertices)
    @test eltype(T) == N
    @test dim(T) == 3
    @test_throws AssertionError Tetrahedron([vertices[1]])
    @test_throws AssertionError Tetrahedron([N[1], N[2], N[3], N[4]])

    vmat = reduce(hcat, vertices)
    Tmat = Tetrahedron(vmat)
    @test Tmat == T

    # constraints_list
    @static if isdefined(@__MODULE__, :Polyhedra)
        c = constraints_list(T)
        @test isequivalent(HPolytope(c), VPolytope(vertices))
    end

    # support vector
    @test σ(ones(3), T) == N[0, 1, 1 / sqrt(2)]

    # membership
    @test zeros(N, 3) ∈ T

    # ispolyhedral
    @test ispolyhedral(T)

    # LazySets#3303
    T = Tetrahedron([N[0, 0, 0], N[0, 1, 0], N[1, 0, 0], N[1, 0, 1]])
    @test zeros(N, 3) ∈ T
end

for N in [Float64, Float32]
    # rand
    @test rand(Tetrahedron; N=N) isa Tetrahedron{N}
end

@test !isoperationtype(Tetrahedron)
