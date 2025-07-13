for N in @tN([Float64, Float32, Rational{Int}])
    # π/2 trigonometric rotation
    b = BallInf(N[1, 2], N(1))
    M = N[0 -1; 1 0]
    # Test Construction
    lm1 = LinearMap(M, b)
    @test lm1.M == M
    @test lm1.X == b
    # Test Dimension
    @test dim(lm1) == 2
    # Test Support Vector
    d = N[1, 1]
    @test σ(d, lm1) == N[-1, 2]
    d = N[-1, 1]
    @test σ(d, lm1) == N[-3, 2]
    d = N[-1, -1]
    @test σ(d, lm1) == N[-3, 0]
    d = N[1, -1]
    @test σ(d, lm1) == N[-1, 0]
    # center
    @test center(lm1) == N[-2, 1]

    # 2D -> 1D Projection
    b = BallInf(N[1, 2], N(1))
    M = N[1 0]
    lm = M * b
    # Test Dimension
    @test dim(lm) == 1
    # Test Support Vector
    d = N[1]
    @test σ(d, lm) == N[2]
    d = N[-1]
    @test σ(d, lm) == N[0]

    # scalar multiplication
    b = BallInf(N[0, 0], N(1))
    lm = N(2) * b
    @test lm == b * N(2)
    # repeated scalar multiplication
    N(2) * lm
    # Test Dimension
    @test dim(lm) == 2
    # Test Support Vector
    d = N[1, 1]
    @test σ(d, lm) == N[2, 2]
    d = N[-1, 1]
    @test σ(d, lm) == N[-2, 2]
    d = N[-1, -1]
    @test σ(d, lm) == N[-2, -2]
    d = N[1, -1]
    @test σ(d, lm) == N[2, -2]

    # Nested construction
    lm1_copy = LinearMap(Matrix{N}(I, 2, 2), lm1)
    @test lm1_copy.M == lm1.M
    @test lm1_copy.X == lm1.X

    # boundedness
    lm1 = ones(N, 2, 2) * Singleton(N[1, 2])  # bounded set
    @test isbounded(lm1) && isboundedtype(typeof(lm1))
    lm2 = zeros(N, 2, 2) * HalfSpace(N[1, 1], N(1))  # zero map
    @test isbounded(lm2) && !isboundedtype(typeof(lm2))
    @test !isbounded(N[2 3; 1 2] * HalfSpace(N[1, 1], N(1)))  # invertible matrix
    @test !isbounded(N[2 3; 0 0] * HalfSpace(N[1, 1], N(1)))  # singular matrix (expensive check)

    # ispolyhedral
    @test ispolyhedral(lm1)
    if N isa AbstractFloat
        lm2 = M * Ball2(N[0, 0], N(1))
        @test !ispolyhedral(lm2)
    end

    # isempty
    @test !isempty(lm)

    # an_element function
    lm = N(2) * BallInf(N[0, 0], N(1))
    an_element(lm)
    @test an_element(lm) ∈ lm

    # check linear map between vector and set
    X = BallInf(N[1], N(1.10445))
    a = N[-1, 2]
    @test a * X isa LinearMap{N,BallInf{N,Vector{N}},N,Matrix{N}}

    # absorbing elements
    for neutral in [ZeroSet{N}(2), EmptySet{N}(2)]
        @test N[0 -1; 1 0] * neutral == neutral
    end

    # vertices_list
    b = BallInf(N[0, 0], N(1))
    M = N[1 2; 3 4]
    vlist = vertices_list(LinearMap(M, b))
    @test ispermutation(vlist, [N[3, 7], N[1, 1], N[-1, -1], N[-3, -7]])
    M = zeros(N, 2, 2)
    vlist = vertices_list(LinearMap(M, b))
    @test vlist == [zeros(N, 2)]

    # test that redundant vertices are removed by default (see #1355)
    X = Interval(N(0), N(1)) × ZeroSet{N}(1)
    M = N[0 1; 0 2]
    @test vertices_list(M * X) == [N[0, 0]]

    # concrete linear map of a LinearMap
    b = BallInf(N[0, 0], N(1))
    M = N[2 3; 1 2]
    L = LinearMap(M, b)
    V = linear_map(M, LinearMap(M, b))

    # concretize
    @test concretize(L) == linear_map(M, b)

    # containment
    L = LineSegment(N[1, 0], N[2, 0])
    # well-conditioned square matrix
    M = N[1 2; 3 4]
    @test N[1, 3] ∈ M * L
    @test N[0, 0] ∉ M * L
    # well-conditioned rectangular matrix
    M = N[1 2; 3 4; 5 6]
    @test N[1, 3, 5] ∈ M * L
    @test N[0, 0, 0] ∉ M * L
    # ill-conditioned square matrix
    M = N[-1 -2; 1 2]
    @test N[-1, 1] ∈ M * L
    @test N[0, 0] ∉ M * L
    # ill-conditioned rectangular matrix
    M = N[-1 -2; 1 2; 5 6]
    @test N[-1, 1, 5] ∈ M * L
    @test N[0, 0, 0] ∉ M * L

    # conversion to polygon in constraint representation
    if N <: AbstractFloat
        Lp = convert(HPolygon, L)
        @test isequivalent(Lp, L)
    end

    # complement
    H = HalfSpace(N[1, 0], N(1)) # x <= 1
    c = complement(N[0 1; 1 0] * H)
    @test c isa UnionSetArray && length(array(c)) == 1
    @test first(array(c)) == HalfSpace(N[0, -1], N(-1)) # complement of y <= 1 is y >= 1
end

for N in [Float64]
    b = BallInf(N[0, 0], N(1))

    if test_suite_polyhedra
        # concrete intersection with lazy linear map
        M = N[2 3; 1 2]
        L = M * b
        Lb = intersection(L, b)
        @test M * an_element(b) ∈ Lb

        # constraints_list
        b = BallInf(N[0, 0], N(1))
        M = N[2 3; 1 2]  # invertible
        lm1 = LinearMap(M, b)
        clist = constraints_list(lm1)
        p1 = HPolygon(clist)
        M = N[2 3; 0 0]  # not invertible
        lm2 = LinearMap(M, b)
        clist = constraints_list(lm2)
        p2 = HPolygon(clist)
        for d in BoxDiagDirections{N}(2)
            @test ρ(d, lm1) ≈ ρ(d, p1)
            @test ρ(d, lm2) ≈ ρ(d, p2)
        end
    end

    # concrete linear map of a LinearMap and membership
    b = BallInf(N[0, 0], N(1))
    M = N[2 3; 1 2]
    L = LinearMap(M, b)
    V = linear_map(M, LinearMap(M, b))
    @test M * M * N[1, 1] ∈ V
end
