using LazySets, Test
using LazySets.ReachabilityBase.Arrays: is_cyclic_permutation
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # constructor from empty vertex list
    P_empty = Polygon{N}()
    @test P_empty isa Polygon{N,Vector{N}} && isempty(P_empty.vertices)

    # constructor from nonempty vertex list
    P = Polygon([N[0, 0], N[0, 2], N[2, 2], N[2, 0], N[1, 1]])

    # convert
    B1 = Hyperrectangle(zeros(N, 2), N[2, 1])
    B2 = Hyperrectangle(zeros(N, 2), N[1, 2])
    B3 = Hyperrectangle(N[0, 3], N[2, 1])
    for U in (UnionSet(B1, B3), UnionSetArray([B1, B3]))
        @test_throws AssertionError convert(Polygon, U)
    end
    for U in (UnionSet(B1, B2), UnionSetArray([B1, B2]))
        P2 = convert(Polygon, U)
        @test is_cyclic_permutation(P2.vertices,
                                    [N[-2, -1], N[-1, -1], N[-1, -2], N[1, -2], N[1, -1], N[2, -1],
                                     N[2, 1], N[1, 1], N[1, 2], N[-1, 2], N[-1, 1], N[-2, 1]])
    end

    # an_element
    @test_throws AssertionError an_element(P_empty)
    x = an_element(P)
    @test x isa Vector{N} && x ∈ P

    # dim
    @test dim(P) == 2

    # boundedness
    @test isbounded(P)

    # isuniversal
    res, w = isuniversal(P, true)
    @test !isuniversal(P) && w ∉ P

    # support vector/function
    d = N[1, 1]
    @test σ(d, P) == N[2, 2]
    @test ρ(d, P) == N(4)

    # isempty
    @test !isempty(P)
    @test isempty(P_empty)

    # convex hull
    @test convex_hull(P) == VPolygon([N[0, 0], N[2, 0], N[2, 2], N[0, 2]])

    # plot_recipe (only check for correct output type)
    x, y = LazySets.plot_recipe(P)
    @test x isa Vector{N} && y isa Vector{N}

    # membership
    P1 = Polygon([N[0, 0], N[2, 1], N[18 // 10, 12 // 10], N[1, 1], N[1 // 2, 3 // 2],
                  N[-1 // 2, 2]])
    P2 = Polygon([N[0, 1], N[2, 1]])
    P3 = Polygon([N[3 // 2, 0], N[3 // 2, 11 // 10]])
    x = N[3 // 2, 1]
    y = N[3 // 2, 3 // 2]
    for P in (P1, P2, P3)
        @test x ∈ P
        @test y ∉ P
    end
    @test x ∈ Polygon([x])
    @test x ∉ Polygon([y])
    @test x ∉ Polygon{N}()
end

# default Float64 constructor
@test Polygon() isa Polygon{Float64,Vector{Float64}}

@test !isoperationtype(Polygon)
@test !isconvextype(Polygon)
@test isboundedtype(Polygon)
