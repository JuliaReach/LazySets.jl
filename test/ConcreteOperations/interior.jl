using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    B = Ball1(zeros(N, 2), N(1 // 10))
    x = N[1 // 10, 1 // 10]
    if N <: AbstractFloat
        @test !is_interior_point(x, B)
    else
        @test_throws ArgumentError is_interior_point(x, B)
        @test !is_interior_point(x, B; ε=N(1 // 100))
    end
    @test !is_interior_point(x, B; ε=N(1))
    @test_throws ArgumentError is_interior_point(x, B; ε=N(0))

    x = N[1 // 10, 1 // 10] .- LazySets._rtol(N)
    if N <: AbstractFloat
        @test !is_interior_point(x, B)
    else
        @test !is_interior_point(x, B; ε=N(1 // 100))
    end
end

for N in @tN([Float64, Float32])
    B = Ball1(zeros(N, 2), N(1 // 10))
    x = N[1 // 10, 1 // 10]
    @test !is_interior_point(x, B; p=N(2))

    x = N[1 // 10, 1 // 10] .- LazySets._rtol(N)
    @test !is_interior_point(x, B; p=N(2))
end
