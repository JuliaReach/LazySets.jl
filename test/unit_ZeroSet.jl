for N in [Float64, Rational{Int}, Float32]
    # random zero set
    rand(ZeroSet)

    Z = ZeroSet{N}(2)
    B = BallInf(ones(N, 2), N(1))

    # support vector
    d = N[1, 0]
    @test σ(d, Z) == zeros(N, 2)

    # testing that the zero set is neutral element for the Minkowski sum
    @test B ⊕ Z == B && Z ⊕ B == B

    cpa = MinkowskiSumArray([B, N(2) * B, N(3) * B])
    @test cpa ⊕ Z == cpa && Z ⊕ cpa == cpa

    # test M-sum of zero set with itself
    @test ZeroSet{N}(2) ⊕ ZeroSet{N}(2) == ZeroSet{N}(2)

    # boundedness
    @test isbounded(Z)

    # element & an_element function
    @test element(Z) ∈ Z
    @test element(Z, 1) == 0
    @test an_element(Z) ∈ Z

    # high and low
    @test high(Z) == zeros(N, 2)
    @test low(Z) == zeros(N, 2)

    # isempty
    @test !isempty(Z)

    # isuniversal
    answer, w = isuniversal(Z, true)
    @test !isuniversal(Z) && !answer && w ∉ Z

    # subset
    z = ZeroSet{N}(1)
    s1 = Singleton(N[0])
    s2 = Singleton(N[2])
    @test z ⊆ s1 && ⊆(z, s1, true)[1]
    subset, point = ⊆(z, s2, true)
    @test z ⊈ s2 && !subset && point ∈ z && point ∉ s2
    @test z ⊆ z && z ⊈ ZeroSet{N}(2)

    # linear map (concrete)
    M = hcat(N[1])
    Mz = linear_map(M, z)
    @test Mz isa ZeroSet && dim(Mz) == 1

    M = N[-1 -2;]
    MZ = linear_map(M, Z)
    @test MZ isa ZeroSet && dim(MZ) == 1

    # translation
    @test translate(Z, N[1, 2]) == Singleton(N[1, 2])

    # base type and element type
    @test basetype(Z) == basetype(typeof(Z)) == ZeroSet
    @test eltype(Z) == eltype(typeof(Z)) == N

    # center
    @test center(Z, 1) == zero(N)

    # rectification
    @test rectify(Z) == Z
end
