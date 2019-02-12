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
#     @test B ⊕ U == U ⊕ B == U ⊕ U == U  # TODO problematic
#     msa = MinkowskiSumArray([B, N(2) * B, N(3) * B])
#     @test msa ⊕ U == U ⊕ msa == U

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

    # an_element
    @test an_element(U) ∈ U

    # norm/radius/diameter functions
    @test_throws ErrorException norm(U)
    @test_throws ErrorException radius(U)
    @test_throws ErrorException diameter(U)
end

# default Float64 constructor
@test Universe(2) == Universe{Float64}(2)
