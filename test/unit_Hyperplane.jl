for N in [Float64, Rational{Int}, Float32]
    # normal constructor
    normal = ones(N, 3)
    hp = Hyperplane(normal, N(5.))

    # dimension
    @test dim(hp) == 3

    # support vector and membership function
    function test_svec(hp, d)
        @test σ(d, hp) ∈ hp
        @test σ(N(2.) * d, hp) ∈ hp
        d2 = N[1., 0., 0.]
        @test_throws ErrorException σ(d2, hp)
        d2 = zeros(N, 3)
        @test σ(d2, hp) ∈ hp
    end
    # tests 1
    normal = ones(N, 3)
    d = ones(N, 3)
    test_svec(Hyperplane(normal, N(5.)), d)
    # tests 2
    normal = zeros(N, 3); normal[3] = N(1.)
    d = zeros(N, 3); d[3] = N(1.)
    test_svec(Hyperplane(normal, N(5.)), d)

    # an_element function and membership function
    @test an_element(hp) ∈ hp
end
