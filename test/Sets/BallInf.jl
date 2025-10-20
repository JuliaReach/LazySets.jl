using LazySets, Test
using LazySets.ReachabilityBase.Arrays: SingleEntryVector, ispermutation
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # invalid inputs
    @test_throws AssertionError BallInf(N[0], N(-1))
    if N <: AbstractFloat
        @test_throws AssertionError BallInf(N[0], N(NaN))
    end

    # 1D BallInf
    b = BallInf(N[0], N(1))
    # Test Dimension
    @test dim(b) == 1
    # Test Support Vector
    d = N[1]
    @test σ(d, b) == N[1]
    d = N[-1]
    @test σ(d, b) == N[-1]

    # 2D BallInf
    b = BallInf(N[0, 0], N(1))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1, 1]
    @test σ(d, b) == N[1, 1]
    d = N[-1, 1]
    @test σ(d, b) == N[-1, 1]
    d = N[-1, -1]
    @test σ(d, b) == N[-1, -1]
    d = N[1, -1]
    @test σ(d, b) == N[1, -1]

    # 2D BallInf not 0-centered
    b = BallInf(N[1, 2], N(1))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1, 1]
    @test σ(d, b) == N[2, 3]
    d = N[-1, 1]
    @test σ(d, b) == N[0, 3]
    d = N[-1, -1]
    @test σ(d, b) == N[0, 1]
    d = N[1, -1]
    @test σ(d, b) == N[2, 1]

    # 2D BallInf radius =/= 1
    b = BallInf(N[0, 0], N(2))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1, 1]
    @test σ(d, b) == N[2, 2]
    d = N[-1, 1]
    @test σ(d, b) == N[-2, 2]
    d = N[-1, -1]
    @test σ(d, b) == N[-2, -2]
    d = N[1, -1]
    @test σ(d, b) == N[2, -2]

    # unicode constructor
    @test □(center(b), radius(b)) == b

    # center
    c = N[1, 2]
    b = BallInf(c, N(2))
    @test center(b) == c && center(b, 1) == N(1) && center(b, 2) == N(2)

    # radius_ball
    @test LazySets.radius_ball(b) == N(2)

    # ball_norm
    @test LazySets.ball_norm(b) == N(Inf)

    # support vector for single entry vector
    svec = σ(SingleEntryVector(2, 3, N(2)), BallInf(zeros(N, 3), N(2)))
    @test -2 <= svec[1] <= 2 && svec[2] == N(2) && -2 <= svec[3] <= 2

    # support function
    B = BallInf(N[1, 2], N(1))
    @test ρ(N[1, 1], B) == N(5)
    @test ρ(N[1, 0], B) == N(2)
    @test ρ(N[0, 1], B) == N(3)
    @test ρ(N[-1, -1], B) == N(-1)
    @test ρ(N[-1, 0], B) == N(0)
    @test ρ(N[-1, 1], B) == N(3)
    # high dimension
    Bh = BallInf(ones(N, 60), N(1))
    @test ρ(ones(N, 60), Bh) == N(120)

    # boundedness
    @test isbounded(b)

    # isoperationtype
    @test !isoperationtype(typeof(b))

    # ispolyhedral
    @test ispolyhedral(b)

    # isempty
    @test !isempty(b)

    # isuniversal
    answer, w = isuniversal(b, true)
    @test !isuniversal(b) && !answer && w ∉ b

    # membership
    b = BallInf(N[1, 1], N(1))
    @test N[0.5, -0.5] ∉ b
    @test N[0.5, 1.5] ∈ b

    # an_element function
    b = BallInf(N[1, 2], N(3))
    @test an_element(b) ∈ b

    # check that vertices_list for zero radius doesn't repeat vertices
    b = BallInf(N[1, 2], N(0))
    vl = vertices_list(b)
    @test vl == [center(b)]

    # low/high/extrema
    b = BallInf(N[1, 2], N(1))
    @test extrema(b) == (low(b), high(b)) == (N[0, 1], N[2, 3])
    @test extrema(b, 1) == (low(b, 1), high(b, 1)) == (N(0), N(2))

    # isflat
    @test !isflat(BallInf(N[1, 1], N(1))) && isflat(BallInf(N[1, 1], N(0)))

    # split
    b = BallInf(N[3, 3], N(1))
    @test split(b, [1, 1]) == [Hyperrectangle(N[3, 3], N[1, 1])]
    @test ispermutation(split(b, [2, 1]),
                        [Hyperrectangle(N[2.5, 3], N[0.5, 1]),
                         Hyperrectangle(N[3.5, 3], N[0.5, 1])])

    # split along the second generator if we see the BallInf as a zonotope
    B2 = BallInf(N[1, 1], N(1))
    Z2 = split(B2, 2)
    @test Z2[1] == Zonotope(N[1, 0.5], N[1 0; 0 0.5])
    @test Z2[2] == Zonotope(N[1, 1.5], N[1 0; 0 0.5])

    # repaeated split method if we see the BallInf as a zonotope
    Z4 = split(B2, [1, 2], [1, 1])
    @test Z4[1] == Zonotope(N[0.5, 0.5], N[0.5 0; 0 0.5])
    @test Z4[2] == Zonotope(N[0.5, 1.5], N[0.5 0; 0 0.5])
    @test Z4[3] == Zonotope(N[1.5, 0.5], N[0.5 0; 0 0.5])
    @test Z4[4] == Zonotope(N[1.5, 1.5], N[0.5 0; 0 0.5])

    # translation
    b = BallInf(N[1, 2], N(1))
    @test translate(b, N[1, 2]) == BallInf(N[2, 4], N(1))
    bb = BallInf(N[0, 0], N(1))
    @test translate!(bb, N[1, 1]) == BallInf(N[1, 1], N(1)) == bb

    # affine map
    M = N[2 1; 0 2]
    B = BallInf(N[0, 0], N(1))
    v = N[1, -1]
    am = affine_map(M, B, v)
    @test ispermutation(vertices_list(am),
                        [N[4, 1], N[0, 1], N[-2, -3], N[2, -3]])
    amv = affine_map(M, B, v; algorithm="vrep") # pass a custom algorithm
    @test amv isa VPolygon && isequivalent(am, amv)

    # area/volume
    B = BallInf(N[0, 0], N(1))
    @test area(B) == volume(B) == N(4)
    B = BallInf(zeros(N, 3), N(2))
    @test area(B) == N(96)
    @test volume(B) == N(64)
    if N <: AbstractFloat
        B = BallInf(zeros(N, 100), N(1 / 2 + 1e-5))
        @test_throws DimensionMismatch area(B)
        @test volume(B) ≈ N(1.0020019812942185)
    end

    # concretize
    B = BallInf(N[0, 0], N(1))
    @test concretize(B) === B

    # constraints iterator
    @test ispermutation(collect(constraints(B)), constraints_list(B))

    # vertices iterator
    @test ispermutation(collect(vertices(B)), vertices_list(B))

    # set difference
    B = BallInf(N[0, 0, 0], N(1))
    @test isempty(difference(B, B))

    # volume
    @test volume(B) ≈ N(8)

    # projection
    b4 = BallInf(N[4, 3, 2, 1], N(2))
    @test project(b4, [2, 4]) == BallInf(N[3, 1], N(2))

    # triangulate (does not work with Float32)
    @static if isdefined(@__MODULE__, :MiniQhull)
        if N != Float32
            B = BallInf(ones(N, 2), N(1))
            Y = triangulate(B; algorithm="delaunay")
            @test Y isa UnionSetArray && length(array(Y)) == 2 &&
                  array(Y)[1] isa VPolytope && array(Y)[2] isa VPolytope
            a = array(Y)
            a1 = [N[0, 0], [2, 0], [2, 2]]
            a2 = [N[0, 0], [0, 2], [2, 2]]
            @test (ispermutation(a[1].vertices, a1) && ispermutation(a[2].vertices, a2)) ||
                  (ispermutation(a[1].vertices, a2) && ispermutation(a[2].vertices, a1))
        end
    end

    # generators
    B = BallInf(N[0, 0], N(1))
    @test ngens(B) == 2
    @test genmat(B) == N[1 0; 0 1]
    @test ispermutation(collect(generators(B)), [N[1, 0], N[0, 1]])
    B_degenerated = BallInf(N[0, 0], N(0))
    @test ngens(B_degenerated) == 0
    gens = genmat(B_degenerated)
    @test gens == Matrix(undef, 2, 0) && gens isa Matrix{N}
    gens = collect(generators(B_degenerated))
    @test isempty(gens) && gens isa Vector{SingleEntryVector{N}}

    # reflect
    @test reflect(b4) == BallInf(N[-4, -3, -2, -1], N(2))

    # scale
    B = BallInf(N[-2, 3], N(1))
    @test scale(N(2), B) == BallInf(N[-4, 6], N(2))
    @test scale(N(-2), B) == BallInf(N[4, -6], N(2))
end

for N in @tN([Float64, Float32])
    # rand
    @test rand(BallInf; N=N) isa BallInf{N}
end

for N in [Float64]
    # concrete Minkowski sum
    b = BallInf(N[1, 2], N(1))
    p = minkowski_sum(b, N[2 0; 0 1] * b)
    @test p isa VPolygon{N} &&
          ispermutation(vertices_list(p), [N[6, 6], N[0, 6], N[0, 2], N[6, 2]])
end
