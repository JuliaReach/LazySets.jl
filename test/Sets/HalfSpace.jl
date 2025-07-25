using LazySets, Test, SparseArrays
using LazySets.ReachabilityBase.Arrays: SingleEntryVector
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # normal constructor
    hs = HalfSpace(ones(N, 3), N(5))

    # corner case: zero normal vector
    @test_throws AssertionError HalfSpace(N[0, 0], N(1))

    # dimension
    @test dim(hs) == 3

    # support function
    hs2 = HalfSpace(N[1, 0], N(1))
    @test ρ(N[2, 0], hs2) == N(2)
    @test ρ(N[-2, 0], hs2) == N(Inf)
    @test ρ(N[1, 1], hs2) == N(Inf)

    # support vector and membership function
    function test_svec(hs, d)
        @test σ(d, hs) ∈ hs
        @test σ(N(2) * d, hs) ∈ hs
        d2 = N[1, 0, 0]
        @test_throws ErrorException σ(d2, hs)
        d2 = zeros(N, 3)
        @test σ(d2, hs) ∈ hs
    end

    # membership with mixed numeric types
    hsm = HalfSpace([-1.0, -1.0], 0.0)
    @test N[1, 1] ∈ hsm
    @test Singleton(N[1, 1]) ⊆ hsm
    @test_throws ArgumentError Singleton(N[1, 1]) ∈ hsm
    @test_throws ArgumentError N[1, 1] ⊆ hsm

    # tests 1
    normal = ones(N, 3)
    d = ones(N, 3)
    test_svec(HalfSpace(normal, N(5)), d)
    # tests 2
    normal = zeros(N, 3)
    normal[3] = N(1)
    d = zeros(N, 3)
    d[3] = 1
    test_svec(HalfSpace(normal, N(5)), d)

    # support vector in other directions throws an error (but see #750)
    # opposite direction
    @test_throws ErrorException σ(N[-1], HalfSpace(N[1], N(1)))
    # any other direction
    @test_throws ErrorException σ(N[1, 1], HalfSpace(N[1, 0], N(1)))

    # boundedness
    @test !isbounded(hs)

    # ispolyhedral
    @test ispolyhedral(hs)

    # universality
    @test !isuniversal(hs)
    res, w = isuniversal(hs, true)
    @test !res && w ∉ hs

    # isempty
    @test !isempty(hs)

    # an_element function and membership function
    @test an_element(hs) ∈ hs

    # constraints list
    @test constraints_list(hs) == [hs]

    # constraints list from matrix-vector representation
    A = N[2 0; 1 3]
    b = N[-1, 1]
    @test constraints_list(A, b) ==
          [HalfSpace(N[2, 0], N(-1)), HalfSpace(N[1, 3], N(1))]

    # constrained dimensions
    @test constrained_dimensions(HalfSpace(N[1, 0, 1], N(1))) == [1, 3]
    @test constrained_dimensions(HalfSpace(N[0, 1, 0], N(1))) == [2]
    # sparse vector
    @test constrained_dimensions(HalfSpace(sparsevec([2], N[1], 3), N(1))) == [2]

    # halfspace_left & halfspace_right
    @test N[1, 2] ∈ halfspace_left(N[1, 1], N[2, 2])
    @test N[2, 1] ∈ halfspace_right(N[1, 1], N[2, 2])

    # translation
    @test translate(hs, N[1, 2, 3]) == HalfSpace(ones(N, 3), N(11))

    # intersection emptiness
    b = BallInf(N[3, 3, 3], N(1))
    empty_intersection, v = is_intersection_empty(b, hs, true)
    @test is_intersection_empty(b, hs) && empty_intersection
    b = BallInf(N[1, 1, 1], N(1))
    empty_intersection, v = is_intersection_empty(b, hs, true)
    @test !is_intersection_empty(b, hs) && !empty_intersection && v ∈ hs
    hs1 = HalfSpace(N[1, 0], N(1)) # x <= 1
    hs2 = HalfSpace(N[-1, 0], N(-2)) # x >= 2
    empty_intersection, v = is_intersection_empty(hs1, hs2, true)
    @test is_intersection_empty(hs1, hs2) && empty_intersection && v == N[]
    hs3 = HalfSpace(N[-1, 0], N(-1)) # x >= 1
    empty_intersection, v = is_intersection_empty(hs1, hs3, true)
    @test !is_intersection_empty(hs1, hs3) && !empty_intersection && v == N[1, 0]
    hs4 = HalfSpace(N[-1, 0], N(0)) # x >= 0
    empty_intersection, v = is_intersection_empty(hs1, hs4, true)
    @test !is_intersection_empty(hs1, hs4) && !empty_intersection && v ∈ hs1 && v ∈ hs4

    # check for tighter constraint
    c1 = HalfSpace(N[1, 0], N(1))
    c2 = HalfSpace(N[2, 0], N(2))
    c3 = HalfSpace(N[2, 0], N(1))
    @test LazySets.is_tighter_same_dir_2D(c1, c2) &&
          LazySets.is_tighter_same_dir_2D(c3, c2) &&
          LazySets.is_tighter_same_dir_2D(c1, c1)
    @test !LazySets.is_tighter_same_dir_2D(c1, c2; strict=true) &&
          LazySets.is_tighter_same_dir_2D(c3, c2; strict=true) &&
          !LazySets.is_tighter_same_dir_2D(c1, c1; strict=true)
    c4 = HalfSpace(N[0, 1], N(2))
    c5 = HalfSpace(N[0, 1], N(1))
    @test !LazySets.is_tighter_same_dir_2D(c4, c5) &&
          LazySets.is_tighter_same_dir_2D(c5, c4)
    @test !LazySets.is_tighter_same_dir_2D(c4, c5; strict=true) &&
          LazySets.is_tighter_same_dir_2D(c5, c4; strict=true)

    # test concrete linear map of a half-space
    H = HalfSpace(N[1, -1], N(0)) # x <= y
    M = N[1 0; 0 0] # non-invertible matrix
    @test_throws ArgumentError linear_map(M, H, algorithm="vrep")
    M = N[2 2; 0 1] # invertible matrix
    @test linear_map(M, H) == HalfSpace(N[0.5, -2.0], N(0.0))
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :CDDLib)
        M = zeros(N, 2, 2) # result is a singleton
        X = linear_map(M, H)
        @test X isa HPolyhedron && isequivalent(X, ZeroSet{N}(2))

        if N isa AbstractFloat
            @test linear_map(N[1 1], H) == Universe{N}(1)
        end
    end

    # projection
    H = HalfSpace(N[1, -1], N(0))  # x <= y
    @test project(H, [1]) == project(H, [2]) == Universe{N}(1)
    @test project(H, [1, 2]) == H

    # conversion of the normal vector
    hs_sev = HalfSpace(SingleEntryVector(2, 3, N(1)), N(1))
    hs_vec = convert(HalfSpace{N,Vector{N}}, hs_sev)
    @test hs_vec.a == N[0, 1, 0] && hs_vec.b == N(1)

    # complement
    @test complement(HalfSpace(N[1, -2], N(3))) == HalfSpace(N[-1, 2], N(-3))

    # complement check
    for (h1, h2, eq) in [(HalfSpace(N[-1], N(0)), HalfSpace(N[1 // 2], N(0)), true),
                         (HalfSpace(N[1, 3], N(1)), HalfSpace(N[-2, -6], N(-2)), true),
                         (HalfSpace(N[1, 3], N(1)), HalfSpace(N[-2, 6], N(-2)), false),
                         (HalfSpace(N[1, 3], N(1)), HalfSpace(N[-2, -6], N(2)), false)]
        @test LazySets.iscomplement(h1, h2) == eq
    end

    # sampling
    for x in sample(H, 10)
        @test x ∈ H
    end

    # permute
    H = HalfSpace(N[1, -2], N(3))
    @test permute(H, 1:2) == H && permute(H, [2, 1]) == HalfSpace(N[-2, 1], N(3))

    # isdisjoint
    H1 = HalfSpace(N[1], N(0))
    H2 = HalfSpace(N[-1], N(-1))
    @test isdisjoint(H1, H2)
    res, w = isdisjoint(H1, H2, true)
    @test res && w isa Vector{N} && w == N[]
    @test !isdisjoint(H1, H1)
    res, w = isdisjoint(H1, H1, true)
    @test !res && w isa Vector{N} && w ∈ H1

    # isfeasible
    clist = [HalfSpace(N[1], N(0))]  # x <= 0
    res, w = isfeasible(clist, true)
    @test isfeasible(clist) && res && w ∈ clist[1]
    clist = [HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-1))]  # x <= 0 && x >= 1
    res, w = isfeasible(clist, true)
    @test !isfeasible(clist) && !res && w isa Vector{N} && isempty(w)

    # remove_redundant_constraints from a list of constraints
    clist = [HalfSpace(N[1], N(1)), HalfSpace(N[1], N(0))]
    clist2 = remove_redundant_constraints(clist)
    res = remove_redundant_constraints!(clist)
    @test res
    @test clist == clist2
    clist = [HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-1)), HalfSpace(N[-1], N(-1))]
    res = remove_redundant_constraints!(clist)
    @test !res
end

for N in @tN([Float64, Float32])
    # rand
    @test rand(HalfSpace; N=N) isa HalfSpace{N}

    # normalization
    hs1 = HalfSpace(N[1e5, 2e5], N(3e5))
    hs2 = normalize(hs1)
    @test norm(hs2.a) ≈ N(1) && hs2.b == hs1.b / norm(hs1.a)
    @test normalize(hs1, N(1)) == HalfSpace(N[1 // 3, 2 // 3], N(1))
    @test normalize(hs1, N(Inf)) == HalfSpace(N[1 // 2, 1], N(3 // 2))

    # distance
    H = HalfSpace(N[1, -1], N(0))  # x <= y
    y = N[1, 1]  # closest point in the half-space
    # point outside
    x = N[2, 0]
    @test distance(x, H) == distance(H, x) ≈ distance(x, y; p=N(2))
    # point at the border
    x = N[1, 1]
    @test distance(x, H) == distance(H, x) ≈ distance(x, y; p=N(2))
    # point strictly inside
    x = N[0, 2]
    @test distance(x, H) == distance(H, x) == N(0)
end

for N in [Float64]
    # rationalization
    H = HalfSpace([1.0, 2.0], 0.0)
    Hr = rationalize(H)
    @test isa(Hr, HalfSpace{Rational{Int},Vector{Rational{Int}}})
    Hr.a == Rational{Int}[1 // 1, 2 // 1]
    Hr = rationalize(BigInt, H)
    @test isa(Hr, HalfSpace{Rational{BigInt},Vector{Rational{BigInt}}})
    Hr.a == Rational{BigInt}[1 // 1, 2 // 1]

    # test robustness of membership function (see LazySets#2312)
    o = N[0.07768723948819561, -0.5762273280928935, 0.28897399484750297, 1.9299362784322858]
    H = HalfSpace(N[-0.09291863543681655, -0.2176689899601838, -0.07453829739226348,
                    0.048948632014371496], N(0.1911363393469332))
    @test o ∈ H

    # tests that require Symbolics
    @static if isdefined(@__MODULE__, :Symbolics)
        # case with only 1 variable
        vars = @variables x
        @test HalfSpace(x <= 2.0, vars) == HalfSpace([1.0], 2.0)

        vars = @variables x y
        @test HalfSpace(2x + 3y < 5) == HalfSpace([2.0, 3.0], 5.0)
        @test HalfSpace(2x + 3y < 5, vars) == HalfSpace([2.0, 3.0], 5.0)
        @test HalfSpace(2x + 3y < 5; N=Int) == HalfSpace([2, 3], 5)

        @test HalfSpace(2x + 3y > 5) == HalfSpace([-2.0, -3.0], -5.0)
        @test HalfSpace(2x + 3y > 5, vars) == HalfSpace([-2.0, -3.0], -5.0)

        @test HalfSpace(2x + 3y ≤ 5) == HalfSpace([2.0, 3.0], 5.0)
        @test HalfSpace(2x + 3y ≤ 5, vars) == HalfSpace([2.0, 3.0], 5.0)

        @test HalfSpace(2x <= 5y - 1) == HalfSpace([2.0, -5.0], -1.0)
        @test HalfSpace(2x ≤ 5y - 1) == HalfSpace([2.0, -5.0], -1.0)

        @test HalfSpace(2x + 3y ≥ 5) == HalfSpace([-2.0, -3.0], -5.0)
        @test HalfSpace(2x + 3y ≥ 5, vars) == HalfSpace([-2.0, -3.0], -5.0)

        @test_throws ArgumentError HalfSpace(2x == 5y, vars)

        # doesn't work because get_vars returns variables [y, x]
        # => both tests below require vars to pass
        @test HalfSpace(2x ≥ 5y - 1, vars) == HalfSpace([-2.0, 5.0], 1.0)
        @test HalfSpace(2x >= 5y - 1, vars) == HalfSpace([-2.0, 5.0], 1.0)

        # test passing a combination of operations
        vars = @variables x[1:2] t
        @test HalfSpace(x[1] <= t, vars) == HalfSpace([1.0, 0.0, -1.0], 0.0)

        # test with sparse variables
        @variables x[1:5]
        @test HalfSpace(2x[1] + 5x[4] <= 10.0, x) == HalfSpace([2.0, 0.0, 0.0, 5.0, 0.0], 10.0)
        @test HalfSpace(2x[1] + 5x[4] >= -10.0 + x[3], x) ==
              HalfSpace([-2.0, 0.0, 1.0, -5.0, 0.0], 10.0)
    end

    # tests that require SymEngine
    @static if isdefined(@__MODULE__, :SymEngine)
        # _ishalfspace
        res = all(LazySets._ishalfspace.([:(x1 <= 0), :(x1 < 0), :(x1 > 0), :(x1 >= 0)]))
        res &= !LazySets._ishalfspace(:(x1 = 0))
        res &= LazySets._ishalfspace(:(2 * x1 <= 4))
        res &= LazySets._ishalfspace(:(6.1 <= 5.3 * f - 0.1 * g))
        res &= !LazySets._ishalfspace(:(2 * x1^2 <= 4))
        res &= !LazySets._ishalfspace(:(x1^2 > 4 * x2 - x3))
        res &= LazySets._ishalfspace(:(x1 > 4 * x2 - x3))
        @test res

        # convert
        H = convert(HalfSpace, :(x1 <= -0.03))
        @test H == HalfSpace([1.0], -0.03)
        H = convert(HalfSpace, :(x1 < -0.03))
        @test H == HalfSpace([1.0], -0.03)
        H = convert(HalfSpace, :(x1 > -0.03))
        @test H == HalfSpace([-1.0], 0.03)
        H = convert(HalfSpace, :(x1 >= -0.03))
        @test H == HalfSpace([-1.0], 0.03)
        H = convert(HalfSpace, :(x1 + x2 <= 2 * x4 + 6))
        @test H == HalfSpace([1.0, 1.0, -2.0], 6.0)
        H = convert(HalfSpace, :(x1 + x2 <= 2 * x4 + 6); vars=SymEngine.Basic[:x1, :x2, :x3, :x4])
        @test H == HalfSpace([1.0, 1.0, 0.0, -2.0], 6.0)
    end
end

# isoperationtype
@test !isoperationtype(HalfSpace)
