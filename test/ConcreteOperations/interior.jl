for N in [Float64, Rational{Int}, Float32]
    P = BallInf(zeros(N, 2), N(1 // 10))
    d = N[1 // 10, 1 // 10]
    if N <: AbstractFloat
        @test !is_interior_point(d, P)
    else
        @test_throws AssertionError is_interior_point(d, P)
        @test !is_interior_point(d, P; ε=1 // 100)
    end
    @test !is_interior_point(d, P; ε=N(1))
    @test_throws AssertionError is_interior_point(d, P; ε=N(0))

    d = N[1 // 10, 1 // 10] .- LazySets._rtol(N)
    if N <: AbstractFloat
        @test is_interior_point(d, P)
    else
        @test !is_interior_point(d, P; ε=1 // 100)
    end
end

# tests that do not work with Rational{Int}
for N in [Float64, Float32]
    P = BallInf(zeros(N, 2), N(1 // 10))
    d = N[1 // 10, 1 // 10]
    @test !is_interior_point(d, P; p=N(2))

    d = N[1 // 10, 1 // 10] .- LazySets._rtol(N)
    @test is_interior_point(d, P; p=N(2))
end
