for N in [Float64, Rational{Int}, Float32]
    # random zonotope
    rand(Zonotope)
    # Test dimension assertion
    @test_throws AssertionError Zonotope(N[1, 1], [N[0 0]; N[1//2 0]; N[0 1//2]; N[0 0]])

    # 1D Zonotope
    z = Zonotope(N[0], Matrix{N}(I, 1, 1))
    # Test Dimension
    @test dim(z) == 1
    # support function
    @test ρ(N[1], z) == ρ(N[-1], z) == N(1)
    # Test Support Vector
    d = N[1]
    @test σ(d, z) == N[1]
    d = N[-1]
    @test σ(d, z) == N[-1]

    # 2D Zonotope
    z = Zonotope(N[0, 0], Matrix{N}(I, 2, 2))
    # Test Dimension
    @test dim(z) == 2
    # support function
    @test ρ(N[1, 0], z) == ρ(N[-1, 0], z) == N(1)
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, z) == N[1, 1] || N[1, -1]
    d = N[-1, 0]
    @test σ(d, z) == N[-1, 1] || N[-1, -1]
    d = N[0, 1]
    @test σ(d, z) == N[1, 1] || N[-1, 1]
    d = N[0, -1]
    @test σ(d, z) == N[1, -1] || N[-1, -1]

    # 2D Zonotope not 0-centered
    z = Zonotope(N[1, 2], Matrix{N}(I, 2, 2))
    # Test Dimension
    @test dim(z) == 2
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, z) == N[2, 3]
    d = N[-1, 0]
    @test σ(d, z) == N[0, 3]
    d = N[0, 1]
    @test σ(d, z) == N[2, 3]
    d = N[0, -1]
    @test σ(d, z) == N[2, 1]

    # zero column in generators
    g = zeros(N, 2, 5)
    g[:, 3] = ones(N, 2)
    g[1, 2] = N(2)
    z = Zonotope(N[1, 2], g)
    @test size(z.generators) == (2, 5)
    zred = remove_zero_generators(z)
    @test size(zred.generators) == (2, 2)

    # boundedness
    @test isbounded(z)

    # is_polyhedral
    @test is_polyhedral(z)

    # isempty
    @test !isempty(z)

    # isuniversal
    answer, w = isuniversal(z, true)
    @test !isuniversal(z) && !answer && w ∉ z

    # an_element function
    @test an_element(z) ∈ z

    # concrete operations
    gens = N[1 1; -1 1]
    Z1 = Zonotope(N[1, 1], gens)
    Z2 = Zonotope(N[-1, 1], Matrix{N}(I, 2, 2))
    A = N[1//2 1; 1 1//2]

    # translation
    @test translate(Z1, N[1, 2]) == Zonotope(N[2, 3], gens)
    Z1c = copy(Z1)
    translate!(Z1c, N[1, 2]) # in-place translation
    @test Z1c == Zonotope(N[2, 3], gens)

    # concrete Minkowski sum
    Z3 = minkowski_sum(Z1, Z2)
    @test Z3.center == N[0, 2]
    @test Z3.generators == N[1 1 1 0; -1 1 0 1]

    # concrete linear map and scale
    Z4 = linear_map(A, Z3)
    @test Z4.center == N[2, 1]
    @test Z4.generators == N[-1//2 3//2 1//2 1; 1//2 3//2 1 1//2]
    Z5 = scale(1//2, Z3)
    @test Z5.center == N[0, 1]
    @test Z5.generators == N[1//2 1//2 1//2 0; -1//2 1//2 0 1//2]

    # in-place linear map
    Zin = convert(Zonotope, BallInf(zeros(N, 2), N(1)))
    Zout = Zonotope(similar(Zin.center), similar(Zin.generators))
    M = N[0 1; -1 0]
    LazySets.linear_map!(Zout, M, Zin)
    @test Zout == Zonotope(N[0, 0], N[0 1; -1 0])

    # in-place scale
    Z5aux = copy(Z3)
    scale!(N(1//2), Z5aux)
    @test isequivalent(Z5, Z5aux)

    # intersection with a hyperplane
    H1 = Hyperplane(N[1, 1], N(3))
    intersection_empty, point = is_intersection_empty(Z1, H1, true)
    @test !is_intersection_empty(Z1, H1) && !intersection_empty
    H2 = Hyperplane(N[1, 1], N(-11))
    @test is_intersection_empty(Z1, H2) && is_intersection_empty(Z1, H2, true)[1]
    @test !is_intersection_empty(H1, Z1)
    @test is_intersection_empty(H2, Z1) && is_intersection_empty(H2, Z1, true)[1]

    # test number of generators
    Z = Zonotope(N[2, 1], N[-1//2 3//2 1//2 1; 1//2 3//2 1 1//2])
    @test ngens(Z) == 4
    @test genmat(Z) == Z.generators
    @test ispermutation(collect(generators(Z)), [genmat(Z)[:, j] for j in 1:ngens(Z)])

    # test order reduction
    Zred1 = reduce_order(Z, 1)
    @test ngens(Zred1) == 2
    @test order(Zred1) == 1
    Zred2 = reduce_order(Z, 2)
    @test ngens(Zred2) == 4
    @test order(Zred2) == 2
    Z = Zonotope(N[2, 1], N[-1//2 3//2 1//2 1 0 1; 1//2 3//2 1 1//2 1 0])
    Zred3 = reduce_order(Z, 2)
    @test ngens(Zred3) == 4
    @test order(Zred3) == 2
    Znogen = Zonotope(N[1, 2], Matrix{N}(undef, 2, 0))
    @test ngens(Znogen) == 0
    @test genmat(Znogen) == Matrix{N}(undef, 2, 0)
    @test collect(generators(Znogen)) == Vector{N}()
    Zs = Zonotope(SVector{2}(Z.center), SMatrix{2, 6}(Z.generators))
    @test reduce_order(Zs, 2) isa Zonotope{N, SVector{2, N}, SMatrix{2, 4, N, 8}}

    # conversion from zonotopic sets
    Z = Zonotope(N[0, 0], hcat(N[1, 1]))
    @test Z == convert(Zonotope, Z) && Z == togrep(Z)
    for AZ in [LineSegment(N[-1, -1], N[1, 1]),
               Hyperrectangle(N[0, 0], N[1, 1]),
               Singleton(N[0, 0]),
               Singleton(sparse(N[0, 0]))]
        Z = Zonotope(center(AZ), genmat(AZ))
        @test convert(Zonotope, AZ) == togrep(AZ) == Z
    end

    # conversion from lazy affine map
    A = N[1 0; 0 1]
    b = N[1, 1]
    B = BallInf(N[0, 0], N(1))
    Z = convert(Zonotope, A * B + b)
    @test Z == Zonotope(N[1, 1], N[1 0; 0 1])

    # test conversion from hyperrectangular sets
    Z = convert(Zonotope, Hyperrectangle(N[2, 3], N[4, 5]))
    @test Z.center == N[2, 3] && diag(Z.generators) == N[4, 5]
    convert(Zonotope, BallInf(N[5, 3], N(2)))
    # flat hyperrectangle
    Z = convert(Zonotope, Hyperrectangle(N[2, 3], N[0, 0]))
    @test Z.center == N[2, 3] && isempty(Z.generators)

    # convert the cartesian product of two hyperrectangles to a zonotope
    h1 = Hyperrectangle(N[1/2],  N[1/2])
    h2 = Hyperrectangle(N[5//2, 9//2],  N[1/2, 1/2])
    H = convert(Hyperrectangle, h1 × h2)
    Z = convert(Zonotope, h1 × h2)
    @test Z ⊆ H && H ⊆ Z

    # same for CartesianProductArray
    Z2 = convert(Zonotope, CartesianProductArray([h1, h2]))
    @test Z == Z2

    # split a zonotope
    Z = Zonotope(N[0, 0], N[1 1; -1 1])
    Z1, Z2 = split(Z, 1) # in this case the splitting is exact
    @test Z1 ⊆ Z && Z2 ⊆ Z
    Z1, Z2, Z3, Z4 = split(Z, [1, 2], [1, 1])
    @test Z1 ⊆ Z && Z2 ⊆ Z && Z3 ⊆ Z && Z4 ⊆ Z
    Z = Zonotope(SVector{2}(N[0, 0]), SMatrix{2, 2}(N[1 1; -1 1]))
    Z1, Z2 = split(Z, 1)
    @test Z1 ⊆ Z && Z2 ⊆ Z

    # converts the cartesian product of two zonotopes to a new zonotope
    Z1 = Zonotope(N[0], hcat(N[1]))
    Z2 = Zonotope(N[1/2], hcat(N[1/2]))
    Z = convert(Zonotope, Z1×Z2)
    @test Z isa Zonotope && Z.center == N[0, 1/2] && Matrix(Z.generators) == N[1 0; 0 1/2]

    # conversion of the lazy linear map of an abstract hyperrectangle to a zonotope
    B = BallInf(N[0], N(1))
    M = hcat(N[1])
    Z = convert(Zonotope, M * B)
    @test Z isa Zonotope && Z.center == N[0] && Z.generators == hcat(N[1])

    # conversion of the lazy linear map of the cartesian product of hyperrectangular
    # sets to a zonotope
    B = BallInf(N[0], N(1))
    M = N[1 0; 0 -1]
    Z = convert(Zonotope, M * (B × B))
    @test Z isa Zonotope && Z.center == N[0, 0] && Z.generators == M

    # same for CPA
    Z2 = convert(Zonotope, M * CartesianProductArray([B, B]))
    @test Z2 == Z

    # list of constraints
    Z = Zonotope(zeros(N, 3), Matrix(N(1)*I, 3, 3))
    B = BallInf(zeros(N, 3), N(1))  # equivalent to Z
    constraints = constraints_list(Z)
    @test constraints isa Vector{<:HalfSpace{N}} && length(constraints) == 6

    # concrete projection returns a zonotope
    πZ12 = project(Z, 1:2)
    @test πZ12 == Zonotope(zeros(N, 2), Matrix(N(1)*I, 2, 2))

    # 1D projection works correctly even with zero generators (#2147)
    Z = convert(Zonotope, BallInf(N[0, 0], N(1)))
    Z2 = project(Z, [1])
    @test Z2 == Zonotope(N[0], hcat(N[1]))

    # removing zero generators in projection is optional
    Z = Zonotope(zeros(N, 3), N[1 0 3; 4 0 6; 0 0 0])
    Z2 = project(Z, [1,2])
    @test Z2 == Zonotope(zeros(N, 2), N[1 3; 4 6])
    Z2 = project(Z, [1,2]; remove_zero_generators=false)
    @test Z2 == Zonotope(zeros(N, 2), N[1 0 3; 4 0 6])

    # remove redundant generators
    zonotopes = [(N[1 2 3 4 5 -1 0;], 1, hcat(16)),
                 (N[1 1 1 1 1 2 -2 0; 0 0 1 1 0 0 0 0; 1 2 0 0 1 2 -2 0], 3, nothing)]
    for (G, nG, G2) in zonotopes
        Z = Zonotope(zeros(N, size(G, 1)), G)
        Z2 = remove_redundant_generators(Z)
        @test ngens(Z2) == nG
        if N<:AbstractFloat || test_suite_polyhedra
            @test isequivalent(Z, Z2)
        end
        if G2 != nothing
            @test genmat(remove_redundant_generators(Z)) == G2
        end

    end
    Z = Zonotope(zeros(N, 3), N[0 0 0; 0 1 0; 0 0 0])
    Z2 = remove_redundant_generators(Z)
    @test center(Z) == center(Z2)
    @test genmat(Z2) == N[0 1 0]'

    # low/high
    Z = Zonotope(N[1, 1], N[1 2 -3; -2 -3 4])
    @test low(Z, 1) == -5 && low(Z, 2) == -8
    @test high(Z, 1) == 7 && high(Z, 2) == 10

    # Minkowski difference (examples taken from Figure 2 in "On computing the
    # Minkowski difference of zonotopes"; the zonotopes Z_s given in the text
    # are not the correct expressions)
    Zm = Zonotope(N[1, 1], N[1 0 1; 0 1 1])
    Zs = Zonotope(N[0, 0], N[1//2 0; -1//8 1//4])
    Zd = Zonotope(N[1, 1], N[1 1//2 0; 1 0 5//8])
    D = minkowski_difference(Zm, Zs)
    @test isequivalent(D, Zd)
    if N != Float32  # not enough precision with Float32
        Zs = Zonotope(N[0, 0], N[1//2 0; -1//2 1//2])
        Zd = Zonotope(N[1, 1], N[1 1//2; 1 0])
        D = minkowski_difference(Zm, Zs)
        @test isequivalent(D, Zd)
    end
    Zs = Zonotope(N[0, 0], N[2 0; -1//2 1//2])
    Zd = EmptySet{N}(2)
    @test isempty(Zd)
end

for N in [Float64]
    # an_element function
    g = zeros(N, 2, 5)
    g[:, 3] = ones(N, 2)
    g[1, 2] = N(2)
    z = Zonotope(N[1, 2], g)
    @test an_element(z) ∈ z

    # conversion to HPolytope
    # 1D
    Z = Zonotope(N[0], Matrix{N}(I, 1, 1))
    P = HPolytope(constraints_list(Z))
    for d in [N[1], N[-1]]
        @test ρ(d, P) == ρ(d, Z)
    end
    # 2D
    Z = Zonotope(N[0, 0], Matrix{N}(I, 2, 2))
    P = HPolytope(constraints_list(Z))
    for d in BoxDiagDirections{N}(2)
        @test ρ(d, P) == ρ(d, Z)
    end
    # sparse matrix (#1468)
    constraints_list(Zonotope(N[0, 0], sparse(N[1 0 ; 0 1])))

    if test_suite_polyhedra
        # constraints_list for generator matrix with a zero row
        Z = Zonotope(N[0, 0], N[2 3; 0 0])
        P = tovrep(HPolygon(constraints_list(Z)))
        @test ispermutation(vertices_list(P), [N[5, 0], [-5, 0]])

        # test that zero generators are ignored (#2147)
        G = spzeros(N, 100, 100)
        G[1, 1] = N(1)
        Z = Zonotope(zeros(N, 100), G)
        v1 = zeros(N, 100)
        v1[1] = N(1)
        v2 = zeros(N, 100)
        v2[1] = N(-1)
        @test ispermutation(vertices_list(Z), [v1, v2])
    end

    # vertices for singleton zonotope (#1881)
    S = Singleton(N[0, 0])
    Z = convert(Zonotope, S)
    @test vertices_list(Z) == [element(S)]

    # test that redundant vertices are removed by default (#1021)
    Z = Zonotope([0., 0.], [1. 0. 1.; 0. 1. 1.])
    vlistZ = vertices_list(Z)
    @test length(vlistZ) == 6
    @test ispermutation(vlistZ, [N[-2, -2], N[0, -2], N[2, 0], N[2, 2], N[0, 2], N[-2, 0]])

    # test 3d zonotope vertex enumeration
    Z = Zonotope([0., 0., 0.], [1. 0. 1.; 0. 1. 1.; 0.1 1. 1.])
    vlistZ = vertices_list(Z)
    @test length(vlistZ) == 8

    # test 2d zonotope generators in positive orthant vertex enumeration
    Z = Zonotope([0., 0.], [1. 0. 1.; 0. 1. 1.])
    vlistZ = vertices_list(Z)
    @test length(vlistZ) == 6
    # redundant vertices are removed automatically
    Z = Zonotope(N[1, 3], N[0 0 1; 0 0 0])
    vlistZ = vertices_list(Z)
    @test length(vlistZ) == 2

    # test 2d zonotope generators in negative orthant vertex enumeration
    Z = Zonotope([0., 0.], -[1. 0. 1.; 0. 1. 1.])
    vlistZ = vertices_list(Z)
    @test length(vlistZ) == 6

    # option to not apply the convex hull operation
    vlistZ = LazySets._vertices_list_iterative(Z.center, Z.generators, apply_convex_hull=false)
    @test length(vlistZ) == 8
    @test ispermutation(convex_hull(vlistZ), [N[-2, -2], N[0, -2], N[2, 0], N[2, 2], N[0, 2], N[-2, 0]])

    # constraints_list correctness check by subset check (requires LP solver)
    B = BallInf(zeros(N, 3), N(1))
    Z = convert(Zonotope, B)
    constraints = constraints_list(Z)
    H = HPolytope(constraints)
    @test H ⊆ B && B ⊆ H

    gens = N[1 1; -1 1]
    Z1 = Zonotope(N[1, 1], gens)
    Z2 = Zonotope(N[-2, -1], Matrix{N}(I, 2, 2))

    # isdisjoint with a hyperplane
    H1 = Hyperplane(N[1, 1], N(3))
    intersection_empty, point = is_intersection_empty(Z1, H1, true)
    @test !intersection_empty && point ∈ Z1 && point ∈ H1
    # zonotope without generators (#2204)
    Z3 = Zonotope(N[0, 0], Matrix{N}(undef, 2, 0))
    @test isdisjoint(Z3, H1)

    # isdisjoint with another zonotope
    result, w = isdisjoint(Z1, Z2, true)
    @test isdisjoint(Z1, Z2) && result && w == N[]
    Z3 = Zonotope(N[2, 1], Matrix{N}(I, 2, 2))
    @test_throws ErrorException isdisjoint(Z1, Z3, true)
    @test !isdisjoint(Z1, Z3)

    # issubset
    Z = Zonotope(N[0, 0], N[1 1; -1 1])
    H1 = Hyperrectangle(low=N[-2, -2], high=N[2, 2])
    H2 = Hyperrectangle(low=N[-2, -2], high=N[2, 0])
    @test issubset(Z, H1)
    @test !issubset(Z, H2)

    # quadratic map
    Z = Zonotope(N[0, 0], N[1 0; 0 1])
    Q1 = N[1/2 0; 0 1/2]
    Q2 = N[0 1/2; 1/2 0]
    # note that there may be repeated generators (though zero generators are removed)
    @test overapproximate(QuadraticMap([Q1, Q2], Z), Zonotope) ==
        Zonotope(N[1//2, 0], N[1//4 1//4 0; 0 0 1])
    Z = Zonotope(N[0, 0], N[1 1; 0 1])
    Q1 = N[1 1; 1 1]
    Q2 = N[-1 0; 0 -1]
    @test overapproximate(QuadraticMap([Q1, Q2], Z), Zonotope) ==
        Zonotope(N[5//2, -3//2], N[1//2 2 4; -1//2 -1 -2])

    # intersection with halfspace
    Z = Zonotope(N[0, 0], N[1 0; 0 1])
    H = HalfSpace(N[1, 0], N(-3//2))
    @test intersection(H, Z) == EmptySet(2)
    H = HalfSpace(N[1, 0], N(3//2))
    @test intersection(H, Z) == Z
    H = HalfSpace(N[1, 0], N(0))
    P = HPolytope([HalfSpace(N[1, 0], N(0)),
                   HalfSpace(N[-1, 0], N(1)),
                   HalfSpace(N[0, 1], N(1)),
                   HalfSpace(N[0, -1], N(1))])
    @test isequivalent(intersection(H, Z), P)
end
