import LazySets.Approximations: decompose,
                                BoxDirections,
                                OctDirections,
                                BoxDiagDirections

for N in [Float64, Rational{Int}, Float32]
    # =====================================
    # Run decompose for different set types
    # =====================================

    function test_directions(set)
        res = σ(N[1, 0], set)[1] == one(N)
        res &= σ(N[0, 1], set)[2] == one(N)
        res &= σ(N[-1, 0], set)[1] == -one(N)
        res &= σ(N[0, -1], set)[2] == -one(N)
        return res
    end
    partition = [[1, 2], [3, 4], [5, 6]]
    b = BallInf(zeros(N, 6), one(N))
    d = decompose(b, partition, HPolygon)
    @test d.array[1] isa HPolygon && test_directions(d.array[1])
    d = decompose(b, partition, Hyperrectangle)
    @test d.array[1] isa Hyperrectangle && test_directions(d.array[1])

    d = decompose(b, partition, LinearMap)
    @test d.array[1] isa LinearMap && test_directions(d.array[1])

    d = decompose(b, [[i] for i in 1:6], Interval)
    @test d.array[1] isa Interval &&
        σ(N[1], d.array[1])[1] == one(N) && σ(N[-1], d.array[1])[1] == -one(N)

    # ===================
    # 1D/3D decomposition
    # ===================

    b1 = Ball1(zeros(N, 7), one(N))
    d = decompose(b1, [[i] for i in 1:7], Hyperrectangle)
    @test length(d.array) == 7
    d = decompose(b1, [[1], 2:3, 4:6, [7]], Hyperrectangle)
    @test length(d.array) == 4

    # ===================
    # template directions
    # ===================

    for dir in [BoxDirections{N}, OctDirections{N}, BoxDiagDirections{N}]
        d = decompose(b, partition, dir)
        @test d isa CartesianProductArray && array(d)[1] isa HPolytope
    end

    # ==========================
    # different options per block
    # ==========================

    block2oa1 = Dict(1 => Hyperrectangle, 2 => HPolygon, 3 => OctDirections{N},
                     4 => Interval)
    block2oa2 = [Hyperrectangle, HPolygon, OctDirections{N}, Interval]
    for block2oa in [block2oa1, block2oa2]  # both Dict and Vector work
        d = decompose(b1, [[1], 2:3, 4:6, [7]], block2oa)
        @test d isa CartesianProductArray && array(d)[1] isa Hyperrectangle &&
              array(d)[2] isa HPolygon && array(d)[3] isa HPolytope &&
              array(d)[4] isa Interval
    end

    # ==================
    # uniform block size
    # ==================

    d = decompose(b, Hyperrectangle; block_size=2)
    @test length(d.array) == 3
    d = decompose(b1, Hyperrectangle; block_size=3)
    @test length(d.array) == 3
end

# tests that do not work with Rational{Int}
for N in [Float64, Float32]
    # =============================
    # Check that Issue #43 is fixed
    # =============================

    CPA = CartesianProductArray
    X = CPA([BallInf(to_N(N, [0.767292, 0.936613]), N(0.1)),
             BallInf(to_N(N, [0.734104, 0.87296]), N(0.1))])
    A = to_N(N, [1.92664 1.00674 1.0731 -0.995149;
        -2.05704 3.48059 0.0317863 1.83481;
        0.990993 -1.97754 0.754192 -0.807085;
        -2.43723 0.782825 -3.99255 3.93324])
    Ω0 = CH(X, A * X)
    dec = decompose(Ω0, [1:2, 3:4], HPolygon)
    dec1 = dec.array[1]

    @test dec1.constraints[1].b ≈ N(2.84042586)
    @test dec1.constraints[2].b ≈ N(4.04708832)
    @test dec1.constraints[3].b ≈ N(-0.667292)
    @test dec1.constraints[4].b ≈ N(-0.836613)

    # =====================
    # ε-close approximation
    # =====================

    function test_directions(set)
        res = σ(N[1, 0], set)[1] == one(N)
        res &= σ(N[0, 1], set)[2] == one(N)
        res &= σ(N[-1, 0], set)[1] == -one(N)
        res &= σ(N[0, -1], set)[2] == -one(N)
        return res
    end
    partition = [[1, 2], [3, 4], [5, 6]]
    b = BallInf(zeros(N, 6), one(N))
    d = decompose(b, partition, HPolygon => N(1e-2))
    @test d.array[1] isa HPolygon && test_directions(d.array[1])
    d = decompose(b, partition, N(1e-2))
    @test d.array[1] isa HPolygon && test_directions(d.array[1])
    d = decompose(b, partition, [N(1e-2), HPolygon => N(1e-2), N(1e-2)])
    @test d.array[1] isa HPolygon && test_directions(d.array[1])
end
