for N in [Float64, Float32, Rational{Int}]
    # -- singletons --
    S = Singleton(N[1, 2])
    @test isequivalent(S, S)
    @test !isequivalent(S, Singleton(N[1, 1]))

    # -- singleton and zonotope --

    S = Singleton(N[1, 2])
    # same center and no generators
    Z = Zonotope(N[1, 2], Matrix{N}(undef, 2, 0))
    @test isequivalent(S, Z) && isequivalent(Z, S)
    # zero generators
    Z = Zonotope(N[1, 2], N[0; 0;;])
    @test isequivalent(S, Z) && isequivalent(Z, S)
    # wrong center
    Z = Zonotope(N[2, 1], Matrix{N}(undef, 2, 0))
    @test !isequivalent(S, Z) && !isequivalent(Z, S)
    # nonzero generators
    Z = Zonotope(N[1, 2], N[1; 0;;])
    @test !isequivalent(S, Z) && !isequivalent(Z, S)
end
