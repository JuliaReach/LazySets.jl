for N in [Float64, Rational{Int}, Float32]
    B = BallInf(ones(N, 2), N(3.))
    H = Hyperrectangle(ones(N, 2), ones(N, 2))
    E = EmptySet{N}()

    # intersection of two sets
    I = Intersection(B, H)

    # dim
    @test dim(I) == 2

    # support vector (currently throws an error)
    @test_throws ErrorException σ(ones(N, 2), I)

    # membership
    @test ∈(ones(N, 2), I) && !∈(N[5., 5.], I)

    # emptiness of intersection
    @test !isempty(I)

    # ---

    # intersection of an array of sets
    IA = IntersectionArray([B, H])

    # dim
    @test dim(IA) == 2

    # support vector (currently throws an error)
    @test_throws ErrorException σ(ones(N, 2), IA)

    # membership
    @test ∈(ones(N, 2), IA) && !∈(N[5., 5.], IA)

    # array getter
    v = Vector{LazySet{N}}(0)
    @test array(IntersectionArray(v)) ≡ v

    # ---

    # absorbing element
    @test I ∩ E == E ∩ I == IA ∩ E == E ∩ IA == E ∩ E == E
end
