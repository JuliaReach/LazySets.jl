for N in [Float64, Rational{Int}, Float32]
    # constructor from empty vertex list
    P = Polygon{N}()
    @test P isa Polygon{N,Vector{N}} && isempty(P.vertices)

    # constructor from nonempty vertex list
    P = Polygon([N[0, 0], N[0, 2], N[2, 2], N[2, 0], N[1, 1]])

    # dim
    @test dim(P) == 2

    # boundedness
    @test isbounded(P)

    # support vector/function
    d = N[1, 1]
    @test σ(d, P) == N[2, 2]
    @test ρ(d, P) == N(4)

    # isempty
    @test !isempty(P)
    @test isempty(Polygon())

    # convex hull
    @test convex_hull(P) == VPolygon([N[0, 0], N[2, 0], N[2, 2], N[0, 2]])

    # plot_recipe (only check for correct output type)
    x, y = LazySets.plot_recipe(P)
    @test x isa Vector{N} && y isa Vector{N}
end

# default Float64 constructor
@test Polygon() isa Polygon{Float64,Vector{Float64}}

@test !isoperationtype(Polygon)
@test !isconvextype(Polygon)
@test isboundedtype(Polygon)
