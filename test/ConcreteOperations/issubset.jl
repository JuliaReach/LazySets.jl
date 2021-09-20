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
end
