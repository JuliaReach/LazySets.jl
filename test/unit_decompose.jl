import LazySets.Approximations.decompose

for N in [Float64, Float32] # TODO Rational{Int}
    # ==============================
    # Check that Issue #43 is fixed
    # ==============================
    const CPA = CartesianProductArray
    X = CPA([BallInf(to_N(N, [0.767292, 0.936613]), to_N(N, 0.1)),
             BallInf(to_N(N, [0.734104, 0.87296]), to_N(N, 0.1))])
    Y = to_N(N, [1.92664 1.00674 1.0731 -0.995149;
        -2.05704 3.48059 0.0317863 1.83481;
        0.990993 -1.97754 0.754192 -0.807085;
        -2.43723 0.782825 -3.99255 3.93324])
    Y = Y * CPA([BallInf(to_N(N, [0.767292, 0.936613]), to_N(N, 0.1)),
                 BallInf(to_N(N, [0.734104, 0.87296]), to_N(N, 0.1))])
    Ω0 = CH(X, Y)
    dec = decompose(Ω0, set_type=HPolygon)
    dec1 = dec.array[1]

    @test dec1.constraints[1].b ≈ to_N(N, 2.84042586)
    @test dec1.constraints[2].b ≈ to_N(N, 4.04708832)
    @test dec1.constraints[3].b ≈ to_N(N, -0.667292)
    @test dec1.constraints[4].b ≈ to_N(N, -0.836613)

    # ======================================
    # Run decompose for different set types
    # ======================================
    function test_directions(set)
        res = σ(N[1, 0], set)[1] == one(N)
        res &= σ(N[0, 1], set)[2] == one(N)
        res &= σ(N[-1, 0], set)[1] == -one(N)
        res &= σ(N[0, -1], set)[2] == -one(N)
        return res
    end
    b = BallInf(zeros(N, 6), one(N))
    d = decompose(b, set_type=HPolygon)
    @test d.array[1] isa HPolygon && test_directions(d.array[1])
    d = decompose(b, set_type=Hyperrectangle)
    @test d.array[1] isa Hyperrectangle && test_directions(d.array[1])
    d = decompose(b, set_type=HPolygon, ɛ=to_N(N, 1e-2))
    @test d.array[1] isa HPolygon && test_directions(d.array[1])

    if N == Float64
    d = decompose(b, set_type=Interval, blocks=ones(Int, 6))
    @test d.array[1] isa Interval &&
        σ(N[1], d.array[1])[1] == one(N) && σ(N[-1], d.array[1])[1] == -one(N)
    end

    # ===================
    # 1D/3D decomposition
    # ===================
    b = Ball2(zeros(N, 7), one(N))
    d = decompose(b, set_type=Hyperrectangle, blocks=ones(Int, 7))
    @test length(d.array) == 7
    d = decompose(b, set_type=Hyperrectangle, blocks=[1,2,3,1])
    @test length(d.array) == 4

    # =======================
    # default block structure
    # =======================
    # even dimension
    b = BallInf(zeros(N, 6), one(N))
    d = decompose(b)
    for ai in array(d)
        @test ai isa Hyperrectangle
    end
    # odd dimension
    b = BallInf(zeros(N, 7), one(N))
    d = decompose(b)
    for ai in array(d)
        @test ai isa Hyperrectangle
    end
    # 1D intervals
    d = decompose(b, set_type=Interval)
    for ai in array(d)
        @test ai isa Interval
    end
end
