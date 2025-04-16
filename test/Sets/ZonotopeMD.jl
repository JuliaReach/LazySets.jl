for N in [Float64, Float32, Rational{Int}]
    # For order 2: 
    #   center ∈ ℝ²,
    #   M is 2×2 and d ∈ ℝ²
    
    # Direct construction via (center, M, d)
    c = [N(1), N(2)]
    M = [N(1) N(0); N(0) N(1)]
    d = [N(1//10), N(2)]
    Zmd = ZonotopeMD(c, M, d)
    @test Zmd isa ZonotopeMD{N}
    @test center(Zmd) == c
    @test genmat(Zmd) == hcat(M, Diagonal(d))
    @test length(Zmd.center) == size(Zmd.M, 1) == length(Zmd.d)
    
    @test length(center(Zmd)) == 2
    @test genmat(Zmd) !== nothing

    # Construction from generator matrix
    G = hcat(M, Diagonal(d))
    Zmd2 = ZonotopeMD(c, G)
    @test Zmd2 == Zmd

    # Generator matrix with wrong shape
    G_wrong = hcat(M, ones(N, 2, 1)) 
    @test_throws AssertionError ZonotopeMD(c, G_wrong)

    # Convert to standard Zonotope 
    Zstd = convert(Zonotope, Zmd)
    @test center(Zstd) == c
    @test genmat(Zstd) == G

    # Cartesian product
    c2 = [N(3), N(4)]
    M2 = [N(2) N(0); N(0) N(2)]
    d2 = [N(3), N(4)]
    Zmd_2 = ZonotopeMD(c2, M2, d2)
    cp_md = cartesian_product(Zmd, Zmd_2)
    cp_std = cartesian_product(Zstd, Zonotope(Zmd_2))
    @test isequivalent(cp_md, cp_std)

    # Apply a linear map
    A = [N(1) N(1); N(0) N(1)]
    Zlm = linear_map(A, Zmd)
    expected_c = A * c
    expected_genmat = A * hcat(M, Diagonal(d))
    @test center(Zlm) == expected_c
    @test genmat(Zlm) == expected_genmat

    # For order 3: 
    #   center ∈ ℝ²,
    #   M is 2×4 and d ∈ ℝ².

    # Direct construction via (center, M, d)
    c = [N(1), N(2)]
    M1 = [N(1) N(0) N(2) N(0); N(0) N(1) N(3) N(0)]
    d = [N(1//10), N(2)]
    Zmd = ZonotopeMD(c, M1, d)
    @test Zmd isa ZonotopeMD{N}
    @test center(Zmd) == c
    @test genmat(Zmd) == hcat(M1, Diagonal(d))
    @test length(Zmd.center) == size(Zmd.M, 1) == length(Zmd.d)
    
    @test length(center(Zmd)) == 2
    @test genmat(Zmd) !== nothing

    # Construction from generator matrix 
    G = hcat(M1, Diagonal(d))
    Zmd2 = ZonotopeMD(c, G)
    @test Zmd2 == Zmd

    G_wrong = hcat(M1, ones(N, 2, 1))  # here p = 3 instead of 4 (2n)
    @test_throws AssertionError ZonotopeMD(c, G_wrong)

    # Convert to standard Zonotope
    Zstd = convert(Zonotope, Zmd)
    @test center(Zstd) == c
    @test genmat(Zstd) == G

    # Cartesian product
    c2 = [N(3), N(4)]
    M2 = [N(2) N(0) N(1) N(0); N(0) N(2) N(3) N(0)]
    d2 = [N(3), N(4)]
    Zmd_2 = ZonotopeMD(c2, M2, d2)
    cp_md = cartesian_product(Zmd, Zmd_2)
    cp_std = cartesian_product(Zstd, Zonotope(Zmd_2))
    @test isequivalent(cp_md, cp_std)

    # Apply a linear map  
    A = [N(1) N(1); N(0) N(1)]
    Zlm = linear_map(A, Zmd)  
    expected_c = A * c
    expected_genmat = A * hcat(M1, Diagonal(d))
    @test center(Zlm) == expected_c
    @test genmat(Zlm) == expected_genmat
end