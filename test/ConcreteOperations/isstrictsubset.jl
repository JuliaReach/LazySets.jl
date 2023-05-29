for N in [Float64, Float32, Rational{Int}]
    H1 = BallInf(N[0, 0], N(1))
    H2 = BallInf(N[0, 0], N(2))

    # strict inclusion
    for (X, Y) in [(H1, H2)]
        res, w = ⊂(X, Y, true)
        @test X ⊂ Y && res && w ∈ Y && w ∉ X
    end

    # no strict inclusion
    for (X, Y) in [(H1, H1), (H2, H1)]
        res, w = ⊂(X, Y, true)
        @test !(X ⊂ Y) && !res && w == N[]
    end
end
