for N in [Float64, Rational{Int}, Float32]
    P = BallInf(zeros(N, 2), N(1//10))
    d = N[1//10, 1//10]
    answer = is_interior_point(d, P)
    @test answer == (N <: Rational)
    @test !is_interior_point(d, P; ε=N(1))
    @test is_interior_point(d, P; ε=N(0))

    d = N[1//10, 1//10] .- LazySets._rtol(N)
    @test is_interior_point(d, P)
end

# tests that do not work with Rational{Int}
for N in [Float64, Float32]
    P = BallInf(zeros(N, 2), N(1//10))
    d = N[1//10, 1//10]
    @test !is_interior_point(d, P; p=N(2))

    d = N[1//10, 1//10] .- LazySets._rtol(N)
    @test is_interior_point(d, P; p=N(2))
end
