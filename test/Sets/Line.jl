using LazySets, Test
using LazySets.ReachabilityBase.Arrays: ispermutation
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # construction
    l1 = Line(; from=N[0, 1], to=N[1, 1]) # two points on the line
    l2 = Line(N[0, 1], N[1, 0]) # point and direction
    @test l1.p ≈ l2.p && l1.d ≈ l2.d  # the lines are the same

    # construction given a 2d direction and offset
    ll = Line(N[0, 1], N(1)) # y == 1
    @test N[0, 1] ∈ ll && N[1, 1] ∈ ll
    ll = Line(N[1, 0], N(1)) # x == 1
    @test N[1, 0] ∈ ll && N[1, 1] ∈ ll
    ll = Line(N[1, -1], N(0)) # x == y
    @test N[0, 0] ∈ ll && N[1, 1] ∈ ll
    @test_throws ArgumentError Line(zeros(N, 2), N(0))

    # direction
    @test direction(l2) == N[1, 0]

    # normalize
    l3 = Line(N[0, 1], N[2, 0])
    @test normalize(l3) == Line(N[0, 1], N[1, 0])
    normalize!(l3)
    @test l3 == Line(N[0, 1], N[1, 0])

    # dimension
    @test dim(l1) == 2

    # isoperationtype
    @test !isoperationtype(Line)

    # support function
    @test ρ(N[0, 1], l1) == N(1)
    @test ρ(N[1, 0], l1) == N(Inf)

    # support vector
    @test σ(N[0, 1], l1) == N[0, 1]
    @test_throws ArgumentError σ(N[1, 0], l1) == N(Inf)

    # boundedness
    @test !isbounded(l1)

    # box_approximation
    @test_throws ArgumentError box_approximation(l1)

    # ispolyhedral
    @test ispolyhedral(l1)

    # universality
    @test !isuniversal(l1)

    # isempty
    @test !isempty(l1)

    # an_element and membership
    an_element(l1) ∈ l1

    # translation
    @test translate(l2, N[0, 1]) == Line(N[0, 2], N[1, 0])

    # distance
    @test distance(N[1, 0], l1) == N(1)
    @test distance(l1, N[1, 0]) == N(1)
    @test distance(Singleton(N[1, 0]), l1) == N(1)
    @test distance(l1, Singleton(N[1, 0])) == N(1)
    @test_throws ArgumentError distance(N[1, 0], l1; p=N(1))

    # concrete linear map special case: singleton result
    l = Line(N[0, 2], N[1, 0])
    @test linear_map(N[0 1; 0 2], l) == Singleton(N[2, 4])

    # projection
    L = Line(N[1, 2, 3], N[1, 0, 0])
    @test project(L, [1]) == Universe{N}(1)
    @test project(L, [2]) == Singleton(N[2])
    @test project(L, [3]) == Singleton(N[3])
    @test project(L, [1, 2]) == Line(N[1, 2], N[1, 0])
    @test project(L, [1, 3]) == Line(N[1, 3], N[1, 0])
    @test project(L, [2, 3]) == Singleton(N[2, 3])
end

for N in @tN([Float64, Float32])
    # rand
    @test rand(Line; N=N) isa Line{N}

    # constraints_list
    l = Line(N[0, 1], N[1, 0])
    @test ispermutation(constraints_list(l),
                        [HalfSpace(N[0, 1], N(1)),    # y <= 1
                         HalfSpace(N[0, -1], N(-1))]) # y >= 1

    # concrete linear map
    mirror = N[-1 0; 0 1]
    l = Line(; from=N[0, 1], to=N[1, 1])
    @test isequivalent(linear_map(mirror, l), l)
    rot = N[0 -1; 1 0] # π/2 ccw rotation
    @test isequivalent(linear_map(rot, l), Line(N[-1, 0], N[0, 1]))

    # construction with normalization
    L = Line(N[0, 1], N[2, 0]; normalize=true)
    @test L == Line(N[0, 1], N[1, 0]) && norm(L.d) ≈ N(1)
    L = Line(; from=N[0, 1], to=N[2, 3], normalize=true)
    @test isequivalent(L, Line(N[0, 1], N[1, 1])) && norm(L.d) ≈ N(1)
    L = Line(N[0, 1], N(1); normalize=true)
    @test isequivalent(L, Line(N[0, 1], N[-1, 0]))
end
