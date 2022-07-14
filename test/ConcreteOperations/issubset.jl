for N in [Float64, Float32, Rational{Int}]

    # interval in union of intervals
    X = Interval(N(0), N(1))

    # included
    YU = Interval(N(1//2), N(2)) ∪ Interval(N(0), N(2//3))
    YU_hr = convert(Hyperrectangle, YU.X) ∪ convert(Hyperrectangle, YU.Y)
    Yarr = UnionSetArray([YU.X, YU.Y])
    for Y in [YU, YU_hr, Yarr]
        res, w = ⊆(X, Y, true)
        @test X ⊆ Y && res && w == N[]
    end

    # not included
    YU = Interval(N(1//2), N(5//6)) ∪ Interval(N(0), N(1//3))
    YU_hr = convert(Hyperrectangle, YU.X) ∪ convert(Hyperrectangle, YU.Y)
    Yarr = UnionSetArray([YU.X, YU.Y])
    for Y in [YU, YU_hr, Yarr]
        res, w = ⊆(X, Y, true)
        @test !(X ⊆ Y) && !res && w ∈ X && w ∉ Y
    end

    # using IA types
    X = IA.interval(N(0), N(1)) # IA
    Y = Interval(N(-1), N(2))
    @test X ⊆ Y
    @test !(Y ⊆ X)

    X = IntervalBox(IA.interval(N(0), N(1)), IA.interval(N(0), N(1)))
    Y = Hyperrectangle(low=[N(-1), N(-1)], high=[N(2), N(2)])
    @test X ⊆ Y
    @test !(Y ⊆ X)

    B = Hyperrectangle(N[0, 0], N[1, 1])
    U = UnionSetArray([Hyperrectangle(N[0, 1], N[3, 1]), Hyperrectangle(N[0, -1], N[3, 1])])
    res, w = ⊆(B, U, true)
    @test B ⊆ U && res && w == N[]
    U = UnionSetArray([Hyperrectangle(N[0, 1], N[3, 1]), Hyperrectangle(N[3, -1], N[3, 1])])
    res, w = ⊆(B, U, true)
    @test !(B ⊆ U) && !res && w ∈ B && w ∉ U
end
