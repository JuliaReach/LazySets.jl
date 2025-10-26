using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # =======
    #  Ball1
    # =======

    # approximation of a centered unit Ball1
    b = Ball1(N[0, 0, 0], N(1))
    @test norm(b) ≈ N(1)  # in the infinity norm (default)
    @test radius(b) ≈ N(1)  # in the infinity norm (default)
    @test diameter(b) ≈ N(2)  # in the infinity norm (default)

    # approximations of a non-centered unit Ball1
    b = Ball1(N[1, 2, 0], N(1))
    @test norm(b) ≈ N(3)  # in the infinity norm (default)
    @test radius(b) ≈ N(1)  # in the infinity norm (default)
    @test diameter(b) ≈ N(2)  # in the infinity norm (default)

    # ================
    #  Hyperrectangle
    # ================

    # metrics in the infinity norm (default)
    h = Hyperrectangle(N[-1, 2], N[0.2, 0.5])
    @test norm(h, Inf) ≈ N(2.5)
    @test radius(h, Inf) ≈ N(0.5)
    @test diameter(h, Inf) ≈ N(1)

    # metrics in the 2-norm
    @test norm(h, 2) ≈ N(sqrt(2.5^2 + 1.2^2))
    @test radius(h, 2) ≈ N(sqrt(0.5^2 + 0.2^2))
    @test diameter(h, 2) ≈ N(2 * sqrt(0.5^2 + 0.2^2))

    # =========
    #  BallInf
    # =========

    # metrics in the infinity norm (default)
    b = BallInf(N[-1, -2, -4], N(0.2))
    @test norm(b, Inf) ≈ N(4.2)
    @test radius(b, Inf) ≈ N(0.2)
    @test diameter(b, Inf) ≈ N(0.4)

    # metrics in the 2-norm
    @test norm(b, 2) ≈ N(sqrt(1.2^2 + 2.2^2 + 4.2^2))
    @test radius(b, 2) ≈ N(sqrt(0.2^2 * 3))
    @test diameter(b, 2) ≈ N(2 * sqrt(0.2^2 * 3))

    # =========
    #  Polygon
    # =========

    # metrics in the infinity norm (default)
    r = N[3 // 10, 2 // 10]
    p = convert(VPolygon, Hyperrectangle(N[0, 1], r))
    @test norm(p, Inf) ≈ N(6 // 5)
    @test radius(p, Inf) ≈ N(3 // 10)
    @test diameter(p, Inf) ≈ N(6 // 10)

    # metrics in the 2-norm
    p33 = VPolytope([N[-1, 0, 0], N[1, 0, 0], N[0, 1, 0]])  # 3D, 3 points
    @test norm(p, 2) ≈ norm(high(p), 2)
    @test norm(p33, 2) == N(1)
    if N <: AbstractFloat
        @test radius(p, 2) ≈ norm(r, 2)
        @test diameter(p, 2) ≈ 2 * norm(r, 2)

        @test radius(p33, 2) == N(1)
        @test diameter(p33, 2) == 2 * N(1)
    end
end
