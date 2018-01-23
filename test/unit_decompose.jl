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
    dec = decompose(Ω0)
    dec1 = dec.sfarray[1]

    @test dec1.constraints_list[1].b ≈ to_N(N, 2.84042586)
    @test dec1.constraints_list[2].b ≈ to_N(N, 4.04708832)
    @test dec1.constraints_list[3].b ≈ to_N(N, -0.667292)
    @test dec1.constraints_list[4].b ≈ to_N(N, -0.836613)

    # ======================================
    # Run decompose for different set types
    # ======================================
    b = BallInf(zeros(N, 6), one(N))
    d = decompose(b, Inf, HPolygon)
    @test d.sfarray[1] isa HPolygon
    d = decompose(b, Inf, Hyperrectangle)
    @test d.sfarray[1] isa Hyperrectangle
    d = decompose(b, to_N(N, 1e-2), Hyperrectangle)
    @test d.sfarray[1] isa HPolygon  # with p not Inf the set type is ignored
    d = decompose(b, to_N(N, 1e-2)) # by default uses HPolygon
    @test d.sfarray[1] isa HPolygon
    d = decompose(b, [to_N(N, 1e-1), to_N(N, 1e-2), to_N(N, 1e-3)])
    @test d.sfarray[1] isa HPolygon
end
