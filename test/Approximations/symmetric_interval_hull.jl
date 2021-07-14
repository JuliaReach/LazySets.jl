for N in [Float64, Rational{Int}, Float32]
    S = Singleton(N[1, -2, 3, -4])
    @test symmetric_interval_hull(S) == Hyperrectangle(zeros(N, 4), N[1, 2, 3, 4])
end
