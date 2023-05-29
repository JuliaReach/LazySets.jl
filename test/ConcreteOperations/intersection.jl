for N in [Float64, Float32, Rational{Int}]
    U = UnionSetArray([Interval(N[-1, 2])])
    I = Interval(N[0, 1])
    @test intersection(U, I) == intersection(I, U) == UnionSetArray([I])
end
