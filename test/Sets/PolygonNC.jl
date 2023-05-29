for N in [Float64, Rational{Int}, Float32]
    P = Polygon([N[0, 0], N[0, 2], N[2, 2], N[2, 0], N[1, 1]])

    # dim
    @test dim(P) == 2

    # boundedness
    @test isbounded(P)

    # support vector/function
    d = N[1, 1]
    @test σ(d, P) == N[2, 2]
    @test ρ(d, P) == N(4)
end

# default Float64 constructor
@test Polygon() isa Polygon{Float64,Vector{Float64}}

@test !isoperationtype(Polygon)
@test !isconvextype(Polygon)
@test isboundedtype(Polygon)
