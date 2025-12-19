using LazySets, Test, LinearAlgebra
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # isoperationtype
    @test !isoperationtype(ZonotopeMD)

    # for order 2:
    #   center ∈ ℝ²,
    #   M is 2×2 and d ∈ ℝ²

    # direct construction via (center, M, d)
    c = [N(1), N(2)]
    M = [N(1) N(0); N(0) N(1)]
    d = [N(1 // 10), N(2)]
    Z = ZonotopeMD(c, M, d)
    @test Z isa ZonotopeMD{N}
    @test center(Z) == c
    G = hcat(M, Diagonal(d))
    @test genmat(Z) == G

    Z1 = ZonotopeMD([N(1)], hcat([N(1)]), [N(2)])  # 1D zonotope

    # dim
    @test dim(Z) == 2

    # ngens
    @test ngens(Z) == 4

    # construction from generator matrix
    Z2 = ZonotopeMD(c, G)
    @test Z2 == Z
    # generator matrix with invalid shape
    G_invalid = hcat(M, ones(N, 2, 1))
    @test_throws AssertionError ZonotopeMD(c, G_invalid)

    # convert to standard zonotope
    Zstd = convert(Zonotope, Z)
    @test center(Zstd) == c
    @test genmat(Zstd) == G

    # high/low
    res = high(Z1)
    @test res isa Vector{N} && res == N[4]
    @test_throws DimensionMismatch high(Z1, 2)
    res = high(Z1, 1)
    @test res isa N && res == N(4)
    res = low(Z1)
    @test res isa Vector{N} && res == N[-2]
    @test_throws DimensionMismatch low(Z1, 2)
    res = low(Z1, 1)
    @test res isa N && res == N(-2)

    # cartesian_product
    c2 = [N(3), N(4)]
    M2 = [N(2) N(0); N(0) N(2)]
    d2 = [N(3), N(4)]
    Z2 = ZonotopeMD(c2, M2, d2)
    cp_md = cartesian_product(Z, Z2)
    @test cp_md isa ZonotopeMD{N}
    cp_std = cartesian_product(Zstd, convert(Zonotope, Z2))
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test isequivalent(cp_md, cp_std)
    end

    # for order 3:
    #   center ∈ ℝ²,
    #   M is 2×4 and d ∈ ℝ².

    # direct construction via (center, M, d)
    c = [N(1), N(2)]
    M = [N(1) N(0) N(2) N(0); N(0) N(1) N(3) N(0)]
    d = [N(1 // 10), N(2)]
    Z = ZonotopeMD(c, M, d)
    @test Z isa ZonotopeMD{N}
    @test center(Z) == c
    @test genmat(Z) == hcat(M, Diagonal(d))

    # dim
    @test dim(Z) == 2

    # ngens
    @test ngens(Z) == 6

    # construction from generator matrix
    G = hcat(M, Diagonal(d))
    Z2 = ZonotopeMD(c, G)
    @test Z2 == Z
    # generator matrix with invalid shape
    G_invalid = hcat(M, ones(N, 2, 1))
    @test_throws AssertionError ZonotopeMD(c, G_invalid)

    # convert to standard zonotope
    Zstd = convert(Zonotope, Z)
    @test center(Zstd) == c
    @test genmat(Zstd) == G

    # cartesian_product
    c2 = [N(3), N(4)]
    M2 = [N(2) N(0) N(1) N(0); N(0) N(2) N(3) N(0)]
    d2 = [N(3), N(4)]
    Z2 = ZonotopeMD(c2, M2, d2)
    cp_md = cartesian_product(Z, Z2)
    @test cp_md isa ZonotopeMD{N}
    cp_std = cartesian_product(Zstd, convert(Zonotope, Z2))
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test isequivalent(cp_md, cp_std)
    end
end

for N in @tN([Float64, Float32])
    # rand
    @test rand(ZonotopeMD, N=N) isa ZonotopeMD{N}
end
