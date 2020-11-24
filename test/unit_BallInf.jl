for N in [Float64, Rational{Int}, Float32]
    # random ball
    rand(BallInf)

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
    d = N[0, -1]
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

    # center
    c = N[1, 2]
    b = BallInf(c, N(2))
    @test center(b) == c && center(b, 1) == N(1) && center(b, 2) == N(2)

    # support vector for single entry vector
    svec = σ(SingleEntryVector(2, 3, N(2)), BallInf(zeros(N, 3), N(2)))
    @test svec[1] ∈ Interval(N[-2, 2]) && svec[2] == N(2) && svec[3] ∈ Interval(N[-2, 2])

    # support function
    B = BallInf(N[1, 2], N(1))
    @test ρ(N[1, 1], B) == N(5)
    @test ρ(N[1, 0], B) == N(2)
    @test ρ(N[0, 1], B) == N(3)
    @test ρ(N[-1, -1], B) == N(-1)
    @test ρ(N[-1, 0], B) == N(0)
    @test ρ(N[-1, 1], B) == N(3)

    # boundedness
    @test isbounded(b)

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

    # high and low
    b = BallInf(N[1, 2], N(1))
    @test high(b) == N[2, 3]
    @test low(b) == N[0, 1]

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

    # affine map
    M = N[2 1; 0 2]
    B = BallInf(N[0, 0], N(1))
    v = N[1, -1]
    am = affine_map(M, B, v)
    @test ispermutation(vertices_list(am),
                        [N[4, 1], N[0, 1], N[-2, -3], N[2, -3]])
    amv = affine_map(M, B, v, algorithm="vrep") # pass a custom algorithm
    @test amv isa VPolygon && isequivalent(am, amv)

    # volume
    B = BallInf(N[0, 0], N(1))
    @test volume(B) ≈ N(4)
    if N <: AbstractFloat
        B = BallInf(zeros(N, 100), N(1/2 + 1e-5))
        @test volume(B) ≈ N(1.0020019812942185)
    end

    # area
    B = BallInf(N[0, 0], N(1))
    @test area(B) ≈ N(4)

    # concretize
    B = BallInf(N[0, 0], N(1))
    @test concretize(B) === B

    # constraints iterator
    @test ispermutation(collect(constraints(B)), constraints_list(B))

    # vertices iterator
    @test ispermutation(collect(vertices(B)), vertices_list(B))

    # conversion to a zonoope
    B = BallInf(N[0, 0], N(1))
    ZB = convert(Zonotope, B)
    @test ZB == Zonotope(N[0, 0], N[1 0; 0 1])
    B = BallInf(N[0, 0], N(0))  # flat case
    ZB = convert(Zonotope, B)
    @test ZB == Zonotope(N[0, 0], Matrix{N}(undef, 2, 0))

   # conversion to a zonotope, static arrays
   B = BallInf(SA[N(0), N(0)], N(1))
   ZB = convert(Zonotope, B)
   @test ZB == Zonotope(SVector{2}(N[0, 0]), SMatrix{2, 2}(N[1 0; 0 1]))
   B = BallInf(SA[N(0), N(0)], N(0))  # flat case
   ZB = convert(Zonotope, B)
   @test ZB == Zonotope(SVector{2}(N[0, 0]), SMatrix{2, 0, N, 0}())
   # specialized method (no prunning)
   B = BallInf(SA[1.0, 2.0], 1.0)
   ZB = LazySets._convert_2D_static(Zonotope, B)
   @test ZB == Zonotope(SA[1.0, 2.0], SA[1.0 0.0; 0.0 1.0])


   # internal function
   B = BallInf(SA[N(0), N(0)], N(1))
   Zs = LazySets._convert_2D_static(Zonotope, B)
   @test Zs == Zonotope(SVector{2}(N[0, 0]), SMatrix{2, 2}(N[1 0; 0 1]))
end

# tests that only work with Float64
for N in [Float64]
    # concrete Minkowski sum
    b = BallInf(N[1, 2], N(1))
    p = minkowski_sum(b, N[2 0; 0 1] * b)
    @test p isa HPolytope{N} && ispermutation(constraints_list(p),
        [HalfSpace(N[0, -1], N(-2)), HalfSpace(N[0, 1], N(6)),
         HalfSpace(N[-1, 0], N(0)), HalfSpace(N[1, 0], N(6))])
end
