using LazySets, Test, LinearAlgebra, SparseArrays
using LazySets.ReachabilityBase.Arrays: ispermutation
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

function isidentical(::Zonotope, ::Zonotope)
    return false
end

function isidentical(Z1::Zonotope{N}, Z2::Zonotope{N}) where {N}
    return Z1.center == Z2.center && Z1.generators == Z2.generators
end

for N in @tN([Float64, Float32, Rational{Int}])
    # auxiliary sets
    # TODO

    # constructor
    c = N[1, 2]
    @test_throws AssertionError Zonotope(c, [N[0 0]; N[1 // 2 0]; N[0 1 // 2]; N[0 0]])
    Z = Zonotope(c, N[1 3; 2 4])
    Z0 = Zonotope(c, N[0 2 1 0 0; 0 0 1 0 0])  # zonotope with zero generators
    Z1 = Zonotope(N[1], hcat(N[1]))  # 1D zonotope
    Z3 = Zonotope(N[1, 1, 1], N[1 1; 2 0; 3 3])  # 3D zonotope
    # constructor from list of generators
    Z2 = Zonotope(c, [N[1, 2], N[3, 4]])
    @test isidentical(Z2, Z)

    # convert
    Z2 = convert(Zonotope, Hyperrectangle(c, N[4, 5]))
    @test isidentical(Z2, Zonotope(c, N[4 0; 0 5]))
    Z = convert(Zonotope, Hyperrectangle(c, N[0, 0]))  # flat hyperrectangle
    @test isidentical(Z2, Zonotope(c, zeros(N, 2, 0)))
    @test_throws AssertionError convert(EmptySet, X)  # TODO define some X



    # an_element
    x = an_element(Z)
    @test x isa Vector{N} && length(x) == 2 && x ∈ Z

    # area
    @test_throws DimensionMismatch area(Z1)
    res = area(Z)
    @test res isa N && res == N()
    res = area(Z3)
    @test res isa N && res == N()

    # chebyshev_center_radius
    @test_throws ArgumentError chebyshev_center_radius(E)

    # complement
    U2 = complement(E)
    @test U2 isa Universe{N} && dim(U2) == 2

    # concretize
    E2 = concretize(E)
    @test isidentical(E, E2)

    # constrained_dimensions
    @test constrained_dimensions(E) == 1:2

    # constraints_list
    @test_throws MethodError constraints_list(E)  # TODO this should maybe change

    # constraints
    @test_throws MethodError constraints(E)  # TODO this should maybe change

    # convex_hull (unary)
    E2 = convex_hull(E)
    @test isidentical(E, E2)

    # copy
    E2 = copy(E)
    @test isidentical(E, E2)

    # diameter
    @test_throws ArgumentError diameter(E, N(1 // 2))
    for res in (diameter(E), diameter(E, Inf), diameter(E, 2))
        @test res isa N && res == N(0)
    end

    # dim
    @test dim(E) == 2

    # eltype
    @test eltype(E) == N
    @test eltype(typeof(E)) == N

    # extrema
    @test_throws ArgumentError extrema(E)
    @test_throws ArgumentError extrema(E, 1)

    # high
    @test_throws ArgumentError high(E)
    @test_throws ArgumentError high(E, 1)

    # isbounded
    @test isbounded(E)

    # isboundedtype
    @test isboundedtype(typeof(E))

    # isconvextype
    @test isconvextype(typeof(E))

    # isempty
    @test isempty(E)
    res, w = isempty(E, true)
    @test res && w isa Vector{N} && isempty(w)

    # isoperation
    @test !isoperation(E)

    # isoperationtype
    @test !isoperationtype(typeof(E))

    # ispolyhedral
    @test !ispolyhedral(E)  # TODO this should maybe change

    # isuniversal
    @test !isuniversal(E)
    res, w = isuniversal(E, true)
    @test !res && w isa Vector{N} && w ∉ E

    # low
    @test_throws ArgumentError low(E)
    @test_throws ArgumentError low(E, 1)

    # norm
    @test_throws ArgumentError norm(E, N(1 // 2))
    for res in (norm(E), norm(E, Inf), norm(E, 2))
        @test res isa N && res == N(0)
    end

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test_throws MethodError polyhedron(E)  # TODO this should maybe change
    end

    # radius
    @test_throws ArgumentError radius(E, N(1 // 2))
    for res in (radius(E), radius(E, Inf), radius(E, 2))
        @test res isa N && res == N(0)
    end

    # rand
    E2 = rand(EmptySet; N=N)
    @test isidentical(E, E2)
    E2 = rand(EmptySet; N=N, dim=3)
    @test isidentical(E3, E2)

    # rectify
    @test rectify(E) == E

    # reflect
    @test reflect(E) == E

    # singleton_list
    res = singleton_list(E)
    T = VERSION < v"1.7" ? Singleton : Singleton{N,Vector{N}}
    @test res isa Vector{T} && isempty(res)

    # tosimplehrep
    @test_throws MethodError tosimplehrep(E)  # TODO this should maybe change

    # triangulate
    @test_throws ArgumentError triangulate(E)

    # triangulate_faces
    @test_throws DimensionMismatch triangulate_faces(E)
    @test_throws ArgumentError triangulate_faces(E3)  # TODO this should maybe change

    # vertices_list
    res = vertices_list(E)
    @test res isa Vector{Vector{N}} && isempty(res)

    # vertices
    res = collect(vertices(E))
    @test res isa Vector{Vector{N}} && isempty(res)

    # volume
    @test volume(E) == N(0)



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

    # ispolyhedral
    @test ispolyhedral(z)

    # isempty
    @test !isempty(z)

    # isuniversal
    answer, w = isuniversal(z, true)
    @test !isuniversal(z) && !answer && w ∉ z

    # membership
    Z = Zonotope(N[1, 2], N[2 1; 1 2])
    @test N[-2, -1] ∈ Z && N[2, 1] ∈ Z && N[4, 5] ∈ Z && N[0, 3] ∈ Z
    @test N[1, 0] ∉ Z && N[3, 2] ∉ Z && N[1, 4] ∉ Z && N[-1, 2] ∉ Z
    # zonotope without generators
    Z = Zonotope(N[1, 2], Matrix{N}(undef, 2, 0))
    @test N[1, 2] ∈ Z && N[1, 1] ∉ Z

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
    Z5 = scale(1 // 2, Z3)
    @test Z5.center == N[0, 1]
    @test Z5.generators == N[1//2 1//2 1//2 0; -1//2 1//2 0 1//2]
    # 1D simplifies to 1 generator
    M = N[-1 1;]
    Z6 = linear_map(M, Z3)
    @test ngens(Z6) == 1 && genmat(Z6) == hcat(N[4])
    # ... unless the map is zero
    M = N[0 0;]
    Z7 = linear_map(M, Z3)
    @test ngens(Z7) == 0 && genmat(Z7) == Matrix{N}(undef, 1, 0)

    # in-place linear map
    Zin = convert(Zonotope, BallInf(zeros(N, 2), N(1)))
    Zout = Zonotope(similar(Zin.center), similar(Zin.generators))
    M = N[0 1; -1 0]
    LazySets.linear_map!(Zout, M, Zin)
    @test Zout == Zonotope(N[0, 0], N[0 1; -1 0])

    # in-place scale
    Z5aux = copy(Z3)
    scale!(N(1 // 2), Z5aux)
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
    for method in [LazySets.ASB10(), LazySets.COMB03(), LazySets.GIR05(), LazySets.SRMB16()]
        Z = Zonotope(N[2, 1], N[-1//2 3//2 1//2 1; 1//2 3//2 1 1//2])
        Zred1 = reduce_order(Z, 1, method)
        @test ngens(Zred1) == 2
        @test order(Zred1) == 1
        Zred2 = reduce_order(Z, 2, method)
        @test ngens(Zred2) == 4
        @test order(Zred2) == 2
        Z = Zonotope(N[2, 1], N[-1//2 3//2 1//2 1 0 1; 1//2 3//2 1 1//2 1 0])
        Zred3 = reduce_order(Z, 2, method)
        @test ngens(Zred3) == 4
        @test order(Zred3) == 2
        Znogen = Zonotope(N[1, 2], Matrix{N}(undef, 2, 0))
        @test ngens(Znogen) == 0
        @test genmat(Znogen) == Matrix{N}(undef, 2, 0)
        @test collect(generators(Znogen)) == Vector{N}()
    end

    @static if isdefined(@__MODULE__, :StaticArrays)
        using StaticArrays: SVector, SMatrix

        # order reduction with static arrays
        Z = Zonotope(N[2, 1], N[-1//2 3//2 1//2 1 0 1; 1//2 3//2 1 1//2 1 0])
        for method in [LazySets.COMB03(), LazySets.GIR05()]
            Zs = Zonotope(SVector{2}(Z.center), SMatrix{2,6}(Z.generators))
            @test reduce_order(Zs, 2, method) isa Zonotope{N,SVector{2,N},SMatrix{2,4,N,8}}
        end
    end

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

    # convert the cartesian product of two hyperrectangles to a zonotope
    h1 = Hyperrectangle(N[1 / 2], N[1 / 2])
    h2 = Hyperrectangle(N[5 // 2, 9 // 2], N[1 / 2, 1 / 2])
    H = convert(Hyperrectangle, h1 × h2)
    Z = convert(Zonotope, h1 × h2)
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test Z ⊆ H && H ⊆ Z
    end

    # same for CartesianProductArray
    Z2 = convert(Zonotope, CartesianProductArray([h1, h2]))
    @test Z == Z2

    # split a zonotope
    Z = Zonotope(N[0, 0], N[1 1; -1 1])
    Z1, Z2 = split(Z, 1) # in this case the splitting is exact
    @test Z1 ⊆ Z && Z2 ⊆ Z
    Z1, Z2, Z3, Z4 = split(Z, [1, 2], [1, 1])
    @test Z1 ⊆ Z && Z2 ⊆ Z && Z3 ⊆ Z && Z4 ⊆ Z
    @static if isdefined(@__MODULE__, :StaticArrays)
        using StaticArrays: SVector, SMatrix

        Z = Zonotope(SVector{2}(N[0, 0]), SMatrix{2,2}(N[1 1; -1 1]))
        Z1, Z2 = split(Z, 1)
        @test Z1 ⊆ Z && Z2 ⊆ Z
    end

    # converts the cartesian product of two zonotopes to a new zonotope
    Z1 = Zonotope(N[0], hcat(N[1]))
    Z2 = Zonotope(N[1 / 2], hcat(N[1 / 2]))
    Z = convert(Zonotope, Z1 × Z2)
    @test Z isa Zonotope && Z.center == N[0, 1 / 2] && Matrix(Z.generators) == N[1 0; 0 1/2]

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

    @static if isdefined(@__MODULE__, :Polyhedra)
        # list of constraints
        Z = Zonotope(zeros(N, 3), Matrix(N(1) * I, 3, 3))
        B = BallInf(zeros(N, 3), N(1))  # equivalent to Z
        clist = constraints_list(Z)
        @test clist isa Vector{<:HalfSpace{N}} && length(clist) == 6
        # test #3209
        Z2 = Zonotope(zeros(N, 3), N[1.0 1 0; 0 1 1; 0 0 0])
        clist = constraints_list(Z2)
        @test clist isa Vector{<:HalfSpace{N}} && length(clist) == 8
    end

    # concrete projection returns a zonotope
    Z = Zonotope(zeros(N, 3), Matrix(N(1) * I, 3, 3))
    πZ12 = project(Z, 1:2)
    @test πZ12 == Zonotope(zeros(N, 2), Matrix(N(1) * I, 2, 2))

    # 1D projection works correctly even with zero generators (#2147)
    Z = convert(Zonotope, BallInf(N[0, 0], N(1)))
    Z2 = project(Z, [1])
    @test Z2 == Zonotope(N[0], hcat(N[1]))

    # removing zero generators in projection is optional
    Z = Zonotope(zeros(N, 3), N[1 0 3; 4 0 6; 0 0 0])
    Z2 = project(Z, [1, 2])
    @test Z2 == Zonotope(zeros(N, 2), N[1 3; 4 6])
    Z2 = project(Z, [1, 2]; remove_zero_generators=false)
    @test Z2 == Zonotope(zeros(N, 2), N[1 0 3; 4 0 6])

    # remove redundant generators
    zonotopes = [(N[1 2 3 4 5 -1 0;], 1, hcat(16)),
                 (N[1 1 1 1 1 2 -2 0; 0 0 1 1 0 0 0 0; 1 2 0 0 1 2 -2 0], 3, nothing)]
    for (G, nG, G2) in zonotopes
        Z = Zonotope(zeros(N, size(G, 1)), G)
        Z2 = remove_redundant_generators(Z)
        @test ngens(Z2) == nG
        if N <: AbstractFloat || isdefined(@__MODULE__, :Polyhedra)
            @test isequivalent(Z, Z2)
        end
        if !isnothing(G2)
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

    # reflect
    Z = Zonotope(N[1, -2], N[2 0; -1//2 1//2])
    @test reflect(Z) == Zonotope(N[-1, 2], N[2 0; -1//2 1//2])

    # permute
    Z = Zonotope(N[-1, -2], N[1 2 3; -4 -5 -6])
    @test permute(Z, [1, 2]) == Z
    @test permute(Z, [2, 1]) == Zonotope(N[-2, -1], N[-4 -5 -6; 1 2 3])

    # norm
    Z = Zonotope(N[1, 2], N[2 1 -2; 1 2 0])
    @test norm(Z, 1) == 11
    @test norm(scale(2.0, Z), 1) == 22.0

    Z2 = Zonotope(N[1, -1], zeros(N, 2, 2))
    @test norm(Z2, 1) == 2

    Z3 = Zonotope(N[1, -1], N[-2 1; 0 -1])
    @test norm(Z3, 1) == 6
end

for N in @tN([Float64, Float32])
    # rand
    @test rand(Zonotope; N=N) isa Zonotope{N}
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
    constraints_list(Zonotope(N[0, 0], sparse(N[1 0; 0 1])))

    @static if isdefined(@__MODULE__, :Polyhedra)
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
    Z = Zonotope([0.0, 0.0], [1.0 0.0 1.0; 0.0 1.0 1.0])
    vlistZ = vertices_list(Z)
    @test length(vlistZ) == 6
    @test ispermutation(vlistZ, [N[-2, -2], N[0, -2], N[2, 0], N[2, 2], N[0, 2], N[-2, 0]])
    c, G = Z.center, Z.generators
    vlist2 = LazySets._vertices_list_2D(c, G; apply_convex_hull=true)
    @test ispermutation(vlistZ, vlist2)
    vlist3 = LazySets._vertices_list_2D(c, G; apply_convex_hull=false)
    @test ispermutation(vlistZ, vlist3)

    @static if isdefined(@__MODULE__, :Polyhedra)
        # test 3d zonotope vertex enumeration
        Z = Zonotope([0.0, 0.0, 0.0], [1.0 0.0 1.0; 0.0 1.0 1.0; 0.1 1.0 1.0])
        vlistZ = vertices_list(Z)
        @test length(vlistZ) == 8
    end

    # test 2d zonotope generators in positive orthant vertex enumeration
    Z = Zonotope([0.0, 0.0], [1.0 0.0 1.0; 0.0 1.0 1.0])
    vlistZ = vertices_list(Z)
    @test length(vlistZ) == 6
    # redundant vertices are removed automatically
    Z = Zonotope(N[1, 3], N[0 0 1; 0 0 0])
    vlistZ = vertices_list(Z)
    @test length(vlistZ) == 2

    # test 2d zonotope generators in negative orthant vertex enumeration
    Z = Zonotope([0.0, 0.0], -[1.0 0.0 1.0; 0.0 1.0 1.0])
    vlistZ = vertices_list(Z)
    @test length(vlistZ) == 6

    # option to not apply the convex hull operation
    vlistZ = LazySets._vertices_list_zonotope_iterative(Z.center, Z.generators;
                                                        apply_convex_hull=false)
    @test length(vlistZ) == 8
    @test ispermutation(convex_hull(vlistZ),
                        [N[-2, -2], N[0, -2], N[2, 0], N[2, 2], N[0, 2], N[-2, 0]])

    # constraints_list correctness check by subset check (requires LP solver)
    B = BallInf(zeros(N, 3), N(1))
    Z = convert(Zonotope, B)
    clist = constraints_list(Z)
    H = HPolytope(clist)
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
    H1 = Hyperrectangle(; low=N[-2, -2], high=N[2, 2])
    H2 = Hyperrectangle(; low=N[-2, -2], high=N[2, 0])
    @test issubset(Z, H1)
    @test !issubset(Z, H2)

    # quadratic map
    Z = Zonotope(N[0, 0], N[1 0; 0 1])
    Q1 = N[1/2 0; 0 1/2]
    Q2 = N[0 1/2; 1/2 0]
    # note that there may be repeated generators (though zero generators are removed)
    @test overapproximate(QuadraticMap([Q1, Q2], Z), Zonotope) ==
          Zonotope(N[1 // 2, 0], N[1//4 1//4 0; 0 0 1])
    Z = Zonotope(N[0, 0], N[1 1; 0 1])
    Q1 = N[1 1; 1 1]
    Q2 = N[-1 0; 0 -1]
    @test overapproximate(QuadraticMap([Q1, Q2], Z), Zonotope) ==
          Zonotope(N[5 // 2, -3 // 2], N[1//2 2 4; -1//2 -1 -2])

    # intersection with halfspace
    Z = Zonotope(N[0, 0], N[1 0; 0 1])
    H = HalfSpace(N[1, 0], N(-3 // 2))
    @test intersection(H, Z) == EmptySet(2)
    H = HalfSpace(N[1, 0], N(3 // 2))
    @test intersection(H, Z) == Z
    H = HalfSpace(N[1, 0], N(0))
    P = HPolytope([HalfSpace(N[1, 0], N(0)),
                   HalfSpace(N[-1, 0], N(1)),
                   HalfSpace(N[0, 1], N(1)),
                   HalfSpace(N[0, -1], N(1))])
    @test isequivalent(intersection(H, Z), P)
end

# isoperationtype
@test !isoperationtype(Zonotope)
