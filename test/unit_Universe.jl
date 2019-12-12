for N in [Float64, Rational{Int}, Float32]
    # random universe
    rand(Universe)

    U = Universe{N}(2)
    B = BallInf(ones(N, 2), N(1))

    # universe is an absorbing element for
    # - convex hull
    @test CH(B, U) == CH(U, B) == CH(U, U) == U
    cha = ConvexHullArray([B, N(2) * B, N(3) * B])
    @test CH(cha, U) == CH(U, cha) == U
    # - union
    @test B ∪ U == U ∪ B == U ∪ U == U
    # - Minkowski sum
    # TODO requires #1099
    # @test B ⊕ U == U ⊕ B == U ⊕ U == U
    # msa = MinkowskiSumArray([B, N(2) * B, N(3) * B])
    # @test msa ⊕ U == U ⊕ msa == U

    # universe is a neutral element for
    # - intersection
    @test B ∩ U == U ∩ B == B
    @test U ∩ U == U
    ia = IntersectionArray([B, N(2) * B, N(3) * B])
    @test ia ∩ U == U ∩ ia == ia

    # dim
    @test dim(U) == 2

    # support function and support vector
    @test ρ(N[0, 1], U) == N(Inf)
    @test ρ(N[0, 0], U) == zero(N)
    @test σ(N[0, 1], U) == N[0, Inf]
    @test σ(N[-1, 2], U) == N[-Inf, Inf]

    # boundedness
    @test !isbounded(U)

    # membership
    @test N[0, 0] ∈ U
    @test_throws AssertionError N[0] ∈ U

    # emptiness
    @test !isempty(U)

    # universality
    @test isuniversal(U) && isuniversal(U, true) == (true, N[])

    # an_element
    @test an_element(U) ∈ U

    # norm/radius/diameter functions
    @test_throws ErrorException norm(U)
    @test_throws ErrorException radius(U)
    @test_throws ErrorException diameter(U)

    # translation
    @test translate(U, N[1, 2]) == U

    # concrete intersection
    @test intersection(B, U) == intersection(U, B) == B
    @test intersection(U, U) == U

    # intersection emptiness
    res, w = isdisjoint(B, U, true)
    @test !isdisjoint(B, U) && !res && w ∈ B && w ∈ U
    res, w = isdisjoint(U, B, true)
    @test !isdisjoint(U, B) && !res && w ∈ B && w ∈ U
    res, w = isdisjoint(U, U, true)
    @test !isdisjoint(U, U) && !res && w ∈ U
    E = EmptySet{N}()
    res, w = isdisjoint(E, U, true)
    @test isdisjoint(E, U) && res && w == N[]
    res, w = isdisjoint(U, E, true)
    @test isdisjoint(U, E) && res && w == N[]

    # subset
    res, w = ⊆(B, U, true)
    @test B ⊆ U && res && w == N[]
    res, w = ⊆(U, B, true)
    @test U ⊈ B && !res && w ∉ B
    res, w = ⊆(U, U, true)
    @test U ⊆ U && res && w == N[]
end

# default Float64 constructor
@test Universe(2) == Universe{Float64}(2)
