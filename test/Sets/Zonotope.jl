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
    Z0 = Zonotope(c, zeros(N, 2, 0))  # zonotope with no generators
    Z1 = Zonotope(N[1], hcat(N[1]))  # 1D zonotope
    Z3 = Zonotope(N[1, 1, 1], N[1 0 0; 0 2 0; 0 0 3])  # 3D zonotope, a hyperrectangle
    # constructor from list of generators
    Z2 = Zonotope(c, [N[1, 2], N[3, 4]])
    @test isidentical(Z2, Z)

    # convert
    # from AbstractZonotope
    Z2 = convert(Zonotope, HParallelotope(N[2 -1; 2 -2], N[1, 2, 1, 1]))
    @test isidentical(Z2, Zonotope(N[-1 // 4, -1 // 2], N[1 -3//4; 1 -3//2]))
    # from AbstractHyperrectangle
    Z2 = convert(Zonotope, Hyperrectangle(c, N[4, 5]))
    @test isidentical(Z2, Zonotope(c, N[4 0; 0 5]))
    Z2 = convert(Zonotope, Hyperrectangle(c, N[0, 0]))  # flat hyperrectangle
    @test isidentical(Z2, Zonotope(c, zeros(N, 2, 0)))
    # TODO static hyperrectangle
    # from CartesianProduct of AbstractHyperrectangles
    # TODO
    # from CartesianProductArray of AbstractHyperrectangles
    # TODO
    # from CartesianProduct of AbstractZonotopes
    # TODO
    # from CartesianProductArray of AbstractZonotopes
    # TODO
    # from LinearMap of AbstractZonotope
    # TODO
    # from LinearMap of CartesianProduct of AbstractHyperrectangles
    # TODO
    # from LinearMap of CartesianProductArray of AbstractHyperrectangles
    # TODO
    # from AbstractAffineMap of AbstractZonotope
    # TODO

    # an_element
    x = an_element(Z)
    @test x isa Vector{N} && x ∈ Z

    # area
    @test_throws DimensionMismatch area(Z1)
    res = area(Z)
    @test res isa N && res == N(8)
    @test_broken area(Z3)  # TODO this should work

    # chebyshev_center_radius
    @static if isdefined(@__MODULE__, :Polyhedra)
        c, r = chebyshev_center_radius(Z)
        @test c isa Vector{N} && c == Z.center
        @test r isa N && n == N(TODO)
    end

    # complement
    if isdefined(@__MODULE__, :Polyhedra) || N <: AbstractFloat
        U = complement(Z3)
        @test U isa UnionSetArray{N} && dim(U) == 3 &&
              ispermutation(array(U),
                            [HalfSpace(N[1, 0, 0], N(0)), HalfSpace(N[-1, 0, 0], N(-2)),
                             HalfSpace(N[0, 1, 0], N(-1)), HalfSpace(N[0, -1, 0], N(-3)),
                             HalfSpace(N[0, 0, 1], N(-2)), HalfSpace(N[0, 0, -1], N(-4))])
    end

    # concretize
    Z2 = concretize(Z)
    @test isidentical(Z, Z2)

    # constrained_dimensions
    @test constrained_dimensions(Z) == 1:2

    if isdefined(@__MODULE__, :Polyhedra) || N <: AbstractFloat
        # constraints_list
        cs = constraints_list(Z3)
        @test ispermutation(cs,
                            [HalfSpace(N[1, 0, 0], N(2)), HalfSpace(N[-1, 0, 0], N(0)),
                             HalfSpace(N[0, 1, 0], N(3)), HalfSpace(N[0, -1, 0], N(1)),
                             HalfSpace(N[0, 0, 1], N(4)), HalfSpace(N[0, 0, -1], N(2))])
        # sparse matrix (#1468)
        Z2 = Zonotope(N[0], sparse(hcat(N[1])))
        res = constraints_list(Z2)
        @test res isa Vector{<:HalfSpace{N}} &&
              ispermutation(res, [HalfSpace(N[1], N(1)), HalfSpace(N[-1], N(1))])

        # constraints
        cs2 = collect(constraints(Z3))
        @test ispermutation(cs2, cs)
    end

    # convex_hull (unary)
    Z2 = convex_hull(Z)
    @test isidentical(Z, Z2)

    # copy
    Z2 = copy(Z)
    @test isidentical(Z, Z2)

    # diameter
    @test_throws ArgumentError diameter(Z, N(1 // 2))
    for res in (diameter(Z), diameter(Z, Inf))
        @test res isa N && res == N(12)
    end
    @test_broken diameter(Z, 2)  # TODO this should work

    # dim
    @test dim(Z) == 2
    @test dim(Z3) == 3

    # eltype
    @test eltype(Z) == N
    @test eltype(typeof(Z)) == N

    # extrema
    res = extrema(Z)
    @test res isa Tuple{Vector{N},Vector{N}} && res[1] == N[-3, -4] && res[2] == N[5, 8]
    @test_throws DimensionMismatch extrema(Z, 3)
    res = extrema(Z, 1)
    @test res isa Tuple{N,N} && res[1] == N(-3) && res[2] == N(5)

    # generators
    @test collect(generators(Z)) == [N[1, 2], N[3, 4]]
    # degenerate case
    gens = collect(generators(Z0))
    @test gens isa AbstractVector{<:AbstractVector{N}} && isempty(gens)

    # genmat
    @test genmat(Z) == Z.generators
    # degenerate case
    gens = genmat(Z0)
    @test gens isa Matrix{N} && gens == Matrix{N}(undef, 2, 0)

    # high
    res = high(Z)
    @test res isa Vector{N} && res == N[5, 8]
    @test_throws DimensionMismatch high(Z, 3)
    res = high(Z, 1)
    @test res isa N && res == N(5)

    # isbounded
    @test isbounded(Z)

    # isboundedtype
    @test isboundedtype(typeof(Z))

    # isconvextype
    @test isconvextype(typeof(Z))

    # isempty
    @test !isempty(Z)
    res, w = isempty(Z, true)
    @test !res && w isa Vector{N} && w ∈ Z

    # isoperation
    @test !isoperation(Z)

    # isoperationtype
    @test !isoperationtype(typeof(Z))

    # ispolyhedral
    @test ispolyhedral(Z)

    # isuniversal
    @test !isuniversal(Z)
    res, w = isuniversal(Z, true)
    @test !res && w isa Vector{N} && w ∉ Z

    # low
    res = low(Z)
    @test res isa Vector{N} && res == N[-3, -4]
    @test_throws DimensionMismatch low(Z, 3)
    res = low(Z, 1)
    @test res isa N && res == N(-3)

    # ngens
    @test ngens(Z) == 2
    # degenerate case
    @test ngens(Z0) == 0

    # norm
    @test_throws ArgumentError norm(Z, N(1 // 2))
    for res in (norm(Z), norm(Z, Inf))
        @test res isa N && res == N(8)
    end
    res = norm(Z, 2)
    @test res isa (N <: AbstractFloat ? N : Float64)
    @test res == norm(N[5, 8], 2)

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        Y = polyhedron(X)
        # TODO
        @test Y isa Polyhedra.Interval{N} && ispermutation(Y.vrep.points.points, [[N(0)], [N(2)]])
    end

    # radius
    @test_throws ArgumentError radius(Z, N(1 // 2))
    for res in (radius(Z), radius(Z, Inf))
        @test res isa N && res == N(6)
    end
    @test_broken radius(Z, 2)  # TODO this should work

    # rectify
    X = rectify(Z)
    @test X isa UnionSetArray{N} && length(X) == 4
    # no precise test here; the `rectify` implementation is tested elsewhere

    # reflect
    Z2 = reflect(Z)
    @test isequivalent(Z2, Zonotope(N[-1, -2], genmat(Z)))

    # singleton_list
    res = singleton_list(Z)
    @test res isa Vector{Singleton{N,Vector{N}}}
    @test ispermutation(res,
                        [Singleton(N[3, 4]), Singleton(N[5, 8]),
                         Singleton(N[-1, 0]), Singleton(N[-3, -4])])

    # tosimplehrep
    A, b = tosimplehrep(Z)
    @test A isa Matrix{N} && b isa Vector{N} && size(A) == (4, 2) && length(b) == 4
    # no precise test here; the `tosimplehrep` implementation is tested elsewhere

    # triangulate
    @static if isdefined(@__MODULE__, :MiniQhull)
        # TODO
        Y = triangulate(Z)
        @test Y isa UnionSetArray{N} && length(Y) == 2
        Y = Y[1]
        @test Y isa VPolytope{N}
        res = vertices_list(Y)
        @test ispermutation(res, [N[0], N[2]])
    end

    # triangulate_faces
    @test_throws DimensionMismatch triangulate_faces(Z)
    @static if isdefined(@__MODULE__, :MiniQhull)
        # TODO
        @test_throws ArgumentError triangulate_faces(Z3)
    end

    # vertices_list
    res = vertices_list(Z)
    @test res isa Vector{Vector{N}} && ispermutation(res, [N[3, 4], N[5, 8], N[-1, 0], N[-3, -4]])
    # degenerate case (#1881)
    res = vertices_list(Z0)
    @test res isa Vector{Vector{N}} && res == [N[1, 2]]
    # zero generators are ignored (#2147)
    Z2 = Zonotope(zeros(N, 2), N[1 0; 0 0])
    @test ispermutation(vertices_list(Z2), [N[1, 0], N[-1, 0]])
    # redundant vertices are removed automatically (#1021)
    Z2 = Zonotope(N[0, 0], N[1 0 1; 0 1 1])
    vlistZ = vertices_list(Z2)
    @test length(vlistZ) == 6
    @test ispermutation(vlistZ, [N[-2, -2], N[0, -2], N[2, 0], N[2, 2], N[0, 2], N[-2, 0]])
    c, G = Z2.center, Z2.generators
    vlist2 = LazySets._vertices_list_2D(c, G; apply_convex_hull=true)
    @test ispermutation(vlistZ, vlist2)
    vlist2 = LazySets._vertices_list_2D(c, G; apply_convex_hull=false)
    @test ispermutation(vlistZ, vlist2)
    # option to not apply the convex-hull operation
    vlistZ = LazySets._vertices_list_zonotope_iterative(Z2.center, Z2.generators;
                                                        apply_convex_hull=false)
    @test length(vlistZ) == 8
    @test ispermutation(convex_hull(vlistZ),
                        [N[-2, -2], N[0, -2], N[2, 0], N[2, 2], N[0, 2], N[-2, 0]])
    # generators in positive resp. negative orthant (i.e., all numbers are positive resp. negative)
    Z2 = Zonotope(N[0, 0], N[1 0 1; 0 1 1])
    for i in 1:2
        if i == 2
            Z2.generators .*= -1
        end
        vlistZ = vertices_list(Z2)
        if i == 1 || N <: AbstractFloat
            @test ispermutation(vlistZ, [N[2, 0], N[2, 2], N[0, 2], N[-2, 0], N[-2, -2], N[0, -2]])
        else
            @test_broken ispermutation(vlistZ,  # TODO this should be fixed
                                       [N[2, 0], N[2, 2], N[0, 2], N[-2, 0], N[-2, -2], N[0, -2]])
        end
    end
    # 3D
    @static if isdefined(@__MODULE__, :Polyhedra)
        vlistZ = vertices_list(Z3)
        @test ispermutation(vlistZ,
                            [N[2, 3, 4], N[2, 3, -2], N[2, -1, 4], N[2, -1, -2],
                             N[0, 3, 4], N[0, 3, -2], N[0, -1, 4], N[0, -1, -2]])
    end

    # vertices
    res = collect(vertices(Z))
    @test res isa Vector{Vector{N}} && ispermutation(res, vertices_list(Z))

    # volume
    @static if isdefined(@__MODULE__, :MiniQhull)
        @test volume(Z) == N(8)  # TODO this should work without Polyhedra
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # TODO revise this loop and merge with above loop

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

    # permute
    Z = Zonotope(N[-1, -2], N[1 2 3; -4 -5 -6])
    @test permute(Z, [1, 2]) == Z
    @test permute(Z, [2, 1]) == Zonotope(N[-2, -1], N[-4 -5 -6; 1 2 3])
end

for N in @tN([Float64, Float32])
    # rand
    Z2 = rand(Zonotope; N=N)
    @test Z2 isa Zonotope{N} && dim(Z2) == 2
    Z2 = rand(Zonotope; N=N, dim=3)
    @test Z2 isa Zonotope{N} && dim(Z2) == 3
    Z2 = rand(Zonotope; N=N, num_generators=5)
    @test Z2 isa Zonotope{N} && dim(Z2) == 2 && ngens(Z2) == 5
end

for N in [Float64]
    # TODO revise this loop

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
