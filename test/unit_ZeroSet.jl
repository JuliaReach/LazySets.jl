for N in [Float64, Rational{Int}, Float32]
    Z = ZeroSet{N}(2)
    B = BallInf(ones(N, 2), N(1.))

    # testing that the zero set is neutral element for the Minkowski sum
    @test B ⊕ Z == B && Z ⊕ B == B

    cpa = MinkowskiSumArray([B, N(2.) * B, N(3.) * B])
    @test cpa ⊕ Z == cpa && Z ⊕ cpa == cpa

    # test M-sum of zero set with itself
    @test ZeroSet{N}(2) ⊕ ZeroSet{N}(2) == ZeroSet{N}(2)

    # an_element function
    @test an_element(Z) ∈ Z

    # subset
    z = ZeroSet{N}(1)
    s1 = Singleton(N[0.])
    s2 = Singleton(N[2.])
    @test ⊆(z, s1) && ⊆(z, s1, true)[1]
    subset, point = ⊆(z, s2, true)
    @test !⊆(z, s2) && !subset && point ∈ z && !(point ∈ s2)
    @test ⊆(z, z) && !⊆(z, ZeroSet{N}(2))
end
