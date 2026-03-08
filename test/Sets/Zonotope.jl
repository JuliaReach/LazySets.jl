using LazySets, Test, LinearAlgebra, SparseArrays
using LazySets.ReachabilityBase.Arrays: ispermutation
using LazySets.ReachabilityBase.Comparison: isapproxzero, set_rtol, _rtol
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
    P = VPolygon([N[3, 4], N[5, 8], N[-1, 0], N[-3, -4]])  # equivalent set
    Xnc = UnionSet(P, BallInf(N[3, 0], N(1)))  # nonconvex set

    # constructor
    c = N[1, 2]
    @test_throws AssertionError Zonotope(c, [N[0 0]; N[1 // 2 0]; N[0 1 // 2]; N[0 0]])
    Z = @inferred Zonotope(c, N[1 3; 2 4])
    Z0 = @inferred Zonotope(c, zeros(N, 2, 0))  # zonotope with no generators
    Z1 = @inferred Zonotope(N[1], hcat(N[1]))  # 1D zonotope
    Z3 = @inferred Zonotope(N[1, 1, 1], N[1 0 0; 0 2 0; 0 0 3])  # 3D zonotope, a hyperrectangle
    # constructor from list of generators
    Z2 = @inferred Zonotope(c, [N[1, 2], N[3, 4]])
    @test isidentical(Z2, Z)

    # convert
    # from AbstractZonotope
    Z2 = @inferred convert(Zonotope, HParallelotope(N[2 -1; 2 -2], N[1, 2, 1, 1]))
    @test isidentical(Z2, Zonotope(N[-1 // 4, -1 // 2], N[1 -3//4; 1 -3//2]))
    # from AbstractHyperrectangle
    r = N[4, 5]
    G = N[4 0; 0 5]
    Z2 = @inferred convert(Zonotope, Hyperrectangle(c, r))
    @test isidentical(Z2, Zonotope(c, G))
    # from AbstractHyperrectangle with flat first dimension
    r2 = N[0, 5]
    Z2 = @inferred convert(Zonotope, Hyperrectangle(c, r2))
    @test isidentical(Z2, Zonotope(c, hcat(r2)))
    # from AbstractHyperrectangle with flat second dimension
    r2 = N[4, 0]
    Z2 = @inferred convert(Zonotope, Hyperrectangle(c, r2))
    @test isidentical(Z2, Zonotope(c, hcat(r2)))
    # from AbstractHyperrectangle flat in both dimensions
    Z2 = @inferred convert(Zonotope, Hyperrectangle(c, N[0, 0]))
    @test isidentical(Z2, Zonotope(c, zeros(N, 2, 0)))
    if isdefined(@__MODULE__, :StaticArrays)
        # from AbstractHyperrectangle with static arrays
        using StaticArrays: SVector, SMatrix

        H = Hyperrectangle(SVector{2}(c), SVector{2}(r))
        Z2a = convert(Zonotope, H)
        @test isidentical(Z2a, Zonotope(SVector{2}(c), SMatrix{2,2}(G)))
        # _convert_static
        Z2 = @inferred LazySets._convert_static(Zonotope, H)
        @test isidentical(Z2, Z2a)
        # flat first dimension
        r2 = N[0, 5]
        Z2 = convert(Zonotope, Hyperrectangle(SVector{2}(c), SVector{2}(r2)))
        @test isidentical(Z2, Zonotope(SVector{2}(c), SMatrix{2,1}(hcat(r2))))
        # flat second dimension
        r2 = N[4, 0]
        Z2 = convert(Zonotope, Hyperrectangle(SVector{2}(c), SVector{2}(r2)))
        @test isidentical(Z2, Zonotope(SVector{2}(c), SMatrix{2,1}(hcat(r2))))
        # flat in both dimensions
        r2 = N[0, 0]
        Z2 = convert(Zonotope, Hyperrectangle(SVector{2}(c), SVector{2}(r2)))
        @test isidentical(Z2, Zonotope(SVector{2}(c), SMatrix{2,0,N,0}()))
    end
    # from AbstractSparsePolynomialZonotope
    Z2 = @inferred convert(Zonotope, convert(SparsePolynomialZonotope, Z))
    @test isidentical(Z2, Z)
    # from ZonotopeMD
    Z2 = @inferred convert(Zonotope, convert(ZonotopeMD, Z))
    @test isequivalent(Z2, Z)
    # from ZonotopeMD with no non-diagonal entry
    Z2 = @inferred convert(Zonotope, ZonotopeMD(N[1, 2], N[1 0; 0 1], N[2, 3]))
    @test isequivalent(Z2, Zonotope(N[1, 2], N[3 0; 0 4]))
    # from CartesianProduct of AbstractHyperrectangles
    H1 = Hyperrectangle([c[1]], [r[1]])
    H2 = Hyperrectangle([c[2]], [r[2]])
    Z2 = @inferred convert(Zonotope, H1 × H2)
    @test isidentical(Z2, Zonotope(c, G))
    # from CartesianProductArray of AbstractHyperrectangles
    Z2 = @inferred convert(Zonotope, CartesianProductArray([H1, H2]))
    @test isidentical(Z2, Zonotope(c, G))
    # from CartesianProduct of AbstractZonotopes
    Z2 = @inferred convert(Zonotope, convert(Zonotope, H1) × convert(Zonotope, H2))
    @test isidentical(Z2, Zonotope(c, G))
    # from CartesianProductArray of AbstractZonotopes
    Z2 = @inferred convert(Zonotope,
                           CartesianProductArray([convert(Zonotope, H1), convert(Zonotope, H2)]))
    @test isidentical(Z2, Zonotope(c, G))
    # from LinearMap of AbstractZonotope
    M = N isa AbstractFloat ? rand(N, 2, 2) : M = N[1 -2; 3 4]
    Z2 = @inferred convert(Zonotope, M * Z)
    @test isidentical(Z2, linear_map(M, Z))
    # from LinearMap of CartesianProduct of AbstractZonotopes
    Z2 = @inferred convert(Zonotope, M * (convert(Zonotope, H1) × convert(Zonotope, H2)))
    @test isidentical(Z2, linear_map(M, Zonotope(c, G)))
    # from LinearMap of CartesianProductArray of AbstractZonotopes
    Z2 = @inferred convert(Zonotope,
                           M *
                           CartesianProductArray([convert(Zonotope, H1), convert(Zonotope, H2)]))
    @test isidentical(Z2, linear_map(M, Zonotope(c, G)))
    # from AbstractAffineMap of AbstractZonotope
    Z2 = @inferred convert(Zonotope, M * Z + c)
    @test isidentical(Z2, affine_map(M, Z, c))

    # an_element
    x = @inferred an_element(Z)
    @test x isa Vector{N} && x ∈ Z

    # area
    @test_throws DimensionMismatch area(Z1)
    if N <: AbstractFloat && VERSION >= v"1.11"
        res = @inferred area(Z)
    else
        res = area(Z)
    end
    @test res isa N && res == N(8)
    if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :GeometryBasics)
        if N <: AbstractFloat && VERSION >= v"1.11"
            res = @inferred area(Z3)
        else
            res = area(Z3)
        end
        @test res == N(88)
        if N isa AbstractFloat
            @test res isa N
        end
    end

    # center
    c2 = @inferred center(Z)
    @test c2 isa AbstractVector{N} && c2 == c
    @test_throws DimensionMismatch center(Z, 3)
    v = @inferred center(Z, 1)
    @test v isa N && v == N(1)

    # complement
    if isdefined(@__MODULE__, :Polyhedra) || N <: AbstractFloat
        if N == Float32 || VERSION < v"1.11"
            X = complement(Z3)
        else
            X = @inferred complement(Z3)
        end
        @test X isa UnionSetArray{N} && dim(X) == 3
        clist = [HalfSpace(N[0, 0, 1], N(-2)), HalfSpace(N[0, 1, 0], N(-1)),
                 HalfSpace(N[1, 0, 0], N(0)), HalfSpace(N[-1, 0, 0], N(-2)),
                 HalfSpace(N[0, -1, 0], N(-3)), HalfSpace(N[0, 0, -1], N(-4))]
        if N <: AbstractFloat
            @test ispermutation(array(X), clist)
        else
            @test all(isequivalent.(array(X), clist))  # rely on order of constraints for this test
        end
    end

    # concretize
    Z2 = @inferred concretize(Z)
    @test isidentical(Z, Z2)

    # constrained_dimensions
    @test (@inferred constrained_dimensions(Z)) == 1:2

    if isdefined(@__MODULE__, :Polyhedra) || N <: AbstractFloat
        # constraints_list
        if N == Float32 || VERSION < v"1.11"
            cs = constraints_list(Z3)
        else
            cs = @inferred constraints_list(Z3)
        end
        clist = [HalfSpace(N[0, 0, -1], N(2)), HalfSpace(N[0, -1, 0], N(1)),
                 HalfSpace(N[-1, 0, 0], N(0)), HalfSpace(N[1, 0, 0], N(2)),
                 HalfSpace(N[0, 1, 0], N(3)), HalfSpace(N[0, 0, 1], N(4))]
        if N <: AbstractFloat
            @test ispermutation(cs, clist)
        else
            @test all(isequivalent.(cs, clist))  # rely on order of constraints for this test
        end
        # sparse matrix (#1468)
        Z2 = Zonotope(N[0], sparse(hcat(N[1])))
        if N == Float32 || VERSION < v"1.11"
            res = constraints_list(Z2)
        else
            res = @inferred constraints_list(Z2)
        end
        @test res isa Vector{<:HalfSpace{N}} &&
              ispermutation(res, [HalfSpace(N[1], N(1)), HalfSpace(N[-1], N(1))])
        # flat zonotope (#3209)
        Z2 = Zonotope(zeros(N, 3), N[1 1 0; 0 1 1; 0 0 0])
        if N == Float32 || VERSION < v"1.11"
            cs2 = constraints_list(Z2)
        else
            cs2 = @inferred constraints_list(Z2)
        end
        if N <: AbstractFloat
            c5 = normalize(HalfSpace(N[1, -1, 0], N(2)))
            c7 = normalize(HalfSpace(N[-1, 1, 0], N(2)))
        else
            c5 = HalfSpace(N[1 // 6, -1 // 6, 0], N(1 // 3))
            c7 = HalfSpace(N[-1 // 6, 1 // 6, 0], N(1 // 3))
        end
        clist = [HalfSpace(N[0, 0, -1], N(0)), HalfSpace(N[0, 0, 1], N(0)),
                 HalfSpace(N[0, -1, 0], N(2)), HalfSpace(N[-1, 0, 0], N(2)), c5,
                 HalfSpace(N[1, 0, 0], N(2)), c7, HalfSpace(N[0, 1, 0], N(2))]
        @test clist isa Vector{<:HalfSpace{N}}
        if N <: AbstractFloat
            @test ispermutation(cs2, clist)
        else
            @test all(isequivalent.(cs2, clist))  # rely on order of constraints for this test
        end

        # constraints
        if N == Float32 || VERSION < v"1.11"
            cs2 = collect(constraints(Z3))
        else
            cs2 = collect(@inferred constraints(Z3))
        end
        @test ispermutation(cs2, cs)
    end

    # convex_hull (unary)
    Z2 = @inferred convex_hull(Z)
    @test isidentical(Z, Z2)

    # copy
    Z2 = @inferred copy(Z)
    @test isidentical(Z, Z2)

    # diameter
    @test_throws ArgumentError diameter(Z, N(1 // 2))
    if N <: AbstractFloat
        for res in ((@inferred diameter(Z)), @inferred diameter(Z, Inf))
            @test res isa N && res == N(12)
        end
        res = @inferred diameter(Z, 2)
        @test res isa N && res == N(14.422205101855956)
    else
        @test_broken @inferred diameter(Z)  # TODO make this type-stable
        for res in (diameter(Z), diameter(Z, Inf))
            @test res isa N && res == N(12)
        end
    end

    # dim
    @test @inferred dim(Z) == 2
    @test @inferred dim(Z3) == 3

    # eltype
    @test (@inferred eltype(Z)) == N
    @test (@inferred eltype(typeof(Z))) == N

    # extrema
    res = @inferred extrema(Z)
    @test res isa Tuple{Vector{N},Vector{N}} && res[1] == N[-3, -4] && res[2] == N[5, 8]
    @test_throws DimensionMismatch extrema(Z, 3)
    res = @inferred extrema(Z, 1)
    @test res isa Tuple{N,N} && res[1] == N(-3) && res[2] == N(5)

    # generators
    @test collect(@inferred generators(Z)) == [N[1, 2], N[3, 4]]
    # degenerate case
    gens = collect(@inferred generators(Z0))
    @test gens isa AbstractVector{<:AbstractVector{N}} && isempty(gens)

    # genmat
    @test (@inferred genmat(Z)) == Z.generators
    # degenerate case
    gens = @inferred genmat(Z0)
    @test gens isa Matrix{N} && gens == Matrix{N}(undef, 2, 0)

    # high
    res = @inferred high(Z)
    @test res isa Vector{N} && res == N[5, 8]
    @test_throws DimensionMismatch high(Z, 3)
    res = @inferred high(Z, 1)
    @test res isa N && res == N(5)

    # isbounded
    @test @inferred isbounded(Z)

    # isboundedtype
    @test @inferred isboundedtype(typeof(Z))

    # isconvex
    @test @inferred isconvex(Z)

    # isconvextype
    @test @inferred isconvextype(typeof(Z))

    # isempty
    @test !(@inferred isempty(Z))
    @test_broken @inferred isempty(Z, true)  # TODO make this type-stable (witness)
    res, w = isempty(Z, true)
    @test !res && w isa Vector{N} && w ∈ Z

    # isoperation
    @test !(@inferred isoperation(Z))

    # isoperationtype
    @test !(@inferred isoperationtype(typeof(Z)))

    # ispolyhedral
    @test @inferred ispolyhedral(Z)

    # ispolyhedraltype
    @test @inferred ispolyhedraltype(typeof(Z))

    # ispolytopic
    @test @inferred ispolytopic(Z)

    # ispolytopictype
    @test @inferred ispolytopictype(typeof(Z))

    # isuniversal
    @test !(@inferred isuniversal(Z))
    @test_broken @inferred isuniversal(Z, true)  # TODO make this type-stable (witness)
    res, w = isuniversal(Z, true)
    @test !res && w isa Vector{N} && w ∉ Z

    # low
    res = @inferred low(Z)
    @test res isa Vector{N} && res == N[-3, -4]
    @test_throws DimensionMismatch low(Z, 3)
    res = @inferred low(Z, 1)
    @test res isa N && res == N(-3)

    # ngens
    @test (@inferred ngens(Z)) == 2
    # degenerate case
    @test (@inferred ngens(Z0)) == 0

    # norm
    @test_throws ArgumentError norm(Z, N(1 // 2))
    if N == Float64 || VERSION < v"1.11"
        res = @inferred norm(Z)
        @test res isa N && res == N(8)
        @static if VERSION >= v"1.11"
            res = @inferred norm(Z, Inf)
        else
            res = norm(Z, Inf)
        end
        @test res isa N && res == N(8)
        @static if VERSION >= v"1.11"
            res = @inferred norm(Z, 2)
        else
            res = norm(Z, 2)
        end
    elseif N == Float32
        for res in ((@inferred norm(Z)), @inferred norm(Z, Inf))
            @test res isa N && res == N(8)
        end
        res = @inferred norm(Z, 2)
    else
        for res in (norm(Z), norm(Z, Inf))
            @test res isa N && res == N(8)
        end
        res = norm(Z, 2)
    end
    @test res isa (N <: AbstractFloat ? N : Float64)
    @test res == norm(N[5, 8], 2)

    # radius
    @test_throws ArgumentError radius(Z, N(1 // 2))
    if N <: AbstractFloat
        for res in ((@inferred radius(Z)), @inferred radius(Z, Inf))
            @test res isa N && res == N(6)
        end
    else
        @test_broken @inferred radius(Z)  # TODO make this type-stable
        for res in (radius(Z), radius(Z, Inf))
            @test res isa N && res == N(6)
        end
    end
    if N <: AbstractFloat
        res = @inferred radius(Z, 2)
        @test res isa N && res == N(7.211102550927978)
    end

    # rectify
    X = rectify(Z)
    @test X isa UnionSetArray{N} && length(X) == 4
    # no precise test here; the `rectify` implementation is tested elsewhere

    # reflect
    Z2 = @inferred reflect(Z)
    @test isequivalent(Z2, Zonotope(N[-1, -2], genmat(Z)))

    # remove_redundant_generators
    Gs = [(N[1 2 3 4 5 -1 0;], hcat(N(16)), false),
          (N[1 1 1 1 1 2 -2 0; 0 0 1 1 0 0 0 0; 1 2 0 0 1 2 -2 0], N[1 6 2; 0 0 2; 2 6 0], false),
          (N[0 0 0; 0 1 0; 0 0 0], N[0 1 0]', false),
          (hcat(N[1, 2]), hcat(N[1, 2]), true)]
    for (G, G2, issame) in Gs
        c2 = zeros(N, size(G, 1))
        Z2a = Zonotope(c2, G)
        Z2b = copy(Z2a)
        Z2c = @inferred remove_redundant_generators(Z2a)
        @test isidentical(Z2a, Z2b)
        if issame
            @test Z2c === Z2a
        else
            if N == Float64
                @test isidentical(Z2c, Zonotope(c2, G2))
            elseif N == Float32
                # imprecision with ±0.0
                @test typeof(Z2c) == typeof(Z2a) && Z2c.center == c2 && Z2c.generators ≈ G2
            else
                # different algorithm results in a different order
                @test typeof(Z2c) == typeof(Z2a) && Z2c.center == c2 &&
                      ispermutation(Vector.(generators(Z2c)), Vector.(eachcol(G2)))
            end
        end
    end

    # remove_zero_generators
    G = zeros(N, 2, 6)
    G[:, 4] = ones(N, 2)
    G[1, 2] = N(2)
    Z2a = Zonotope(c, G)
    Z2b = copy(Z2a)
    Z2c = @inferred remove_zero_generators(Z2a)
    @test isidentical(Z2a, Z2b)
    @test isidentical(Z2c, Zonotope(c, N[2 1; 0 1]))

    # singleton_list
    @static if VERSION >= v"1.11"
        res = @inferred singleton_list(Z)
    else
        res = singleton_list(Z)
    end
    @test res isa Vector{Singleton{N,Vector{N}}}
    @test ispermutation(res,
                        [Singleton(N[3, 4]), Singleton(N[5, 8]),
                         Singleton(N[-1, 0]), Singleton(N[-3, -4])])

    # togrep
    Z2 = @inferred togrep(Z)
    @test isidentical(Z2, Z)

    # tosimplehrep
    if (N == Float32) || (VERSION < v"1.11")
        A, b = tosimplehrep(Z)
    else
        A, b = @inferred tosimplehrep(Z)
    end
    @test A isa Matrix{N} && b isa Vector{N} && size(A) == (4, 2) && length(b) == 4
    # no precise test here; the `tosimplehrep` implementation is tested elsewhere

    # triangulate_faces
    @test_throws DimensionMismatch triangulate_faces(Z)
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :GeometryBasics)
        res = @inferred triangulate_faces(Z3)
        @test res isa Tuple{Matrix{Float32},Vector{Tuple{Int,Int,Int}}}
        @test length(res[2]) == 12
        # no precise test here; the `triangulate_faces` implementation is tested elsewhere
    end

    # vertices_list
    @static if VERSION >= v"1.11"
        res = @inferred vertices_list(Z)
    else
        res = vertices_list(Z)
    end
    @test res isa Vector{Vector{N}} && ispermutation(res, [N[3, 4], N[5, 8], N[-1, 0], N[-3, -4]])
    # degenerate case (#1881)
    @static if VERSION >= v"1.11"
        res = @inferred vertices_list(Z0)
    else
        res = vertices_list(Z0)
    end
    @test res isa Vector{Vector{N}} && res == [N[1, 2]]
    # zero generators are ignored (#2147)
    Z2 = Zonotope(zeros(N, 2), N[1 0; 0 0])
    @static if VERSION >= v"1.11"
        res = @inferred vertices_list(Z2)
    else
        res = vertices_list(Z2)
    end
    @test ispermutation(res, [N[1, 0], N[-1, 0]])
    # redundant vertices are removed automatically (#1021)
    Z2 = Zonotope(N[0, 0], N[1 0 1; 0 1 1])
    @static if VERSION >= v"1.11"
        vlistZ = @inferred vertices_list(Z2)
    else
        vlistZ = vertices_list(Z2)
    end
    @test length(vlistZ) == 6
    @test ispermutation(vlistZ, [N[-2, -2], N[0, -2], N[2, 0], N[2, 2], N[0, 2], N[-2, 0]])
    c2, G = Z2.center, Z2.generators
    @test_broken @inferred LazySets._vertices_list_2D(c2, G; apply_convex_hull=true)  # TODO make this type-stable
    vlist2 = LazySets._vertices_list_2D(c2, G; apply_convex_hull=true)
    @test ispermutation(vlistZ, vlist2)
    vlist2 = LazySets._vertices_list_2D(c2, G; apply_convex_hull=false)
    @test ispermutation(vlistZ, vlist2)
    # option to not apply the convex-hull operation
    vlistZ = @inferred LazySets._vertices_list_zonotope_iterative(Z2.center, Z2.generators;
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
        @static if VERSION >= v"1.11"
            vlistZ = @inferred vertices_list(Z2)
        else
            vlistZ = vertices_list(Z2)
        end
        @test ispermutation(vlistZ, [N[2, 0], N[2, 2], N[0, 2], N[-2, 0], N[-2, -2], N[0, -2]])
    end
    # 3D
    @static if isdefined(@__MODULE__, :Polyhedra)
        @static if VERSION >= v"1.11"
            vlistZ = @inferred vertices_list(Z3)
        else
            vlistZ = vertices_list(Z3)
        end
        @test ispermutation(vlistZ,
                            [N[2, 3, 4], N[2, 3, -2], N[2, -1, 4], N[2, -1, -2],
                             N[0, 3, 4], N[0, 3, -2], N[0, -1, 4], N[0, -1, -2]])
    end

    # vertices
    @static if VERSION >= v"1.11"
        res = collect(@inferred vertices(Z))
    else
        res = collect(vertices(Z))
    end
    @test res isa Vector{Vector{N}} && ispermutation(res, vertices_list(Z))

    # volume (part 1)
    if N == Float64 && VERSION >= v"1.11"
        res = @inferred volume(Z)
    else
        res = volume(Z)
    end
    @test res isa N && res == N(8)

    # affine_map
    @test_throws DimensionMismatch affine_map(ones(N, 2, 1), Z, N[1, 1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 2), Z, N[1])
    Z2 = @inferred affine_map(ones(N, 1, 2), Z, N[1])
    @test isidentical(Z2, Zonotope(N[4], hcat(N(10))))
    Z2 = @inferred affine_map(ones(N, 3, 2), Z, N[1, 2, 3])
    @test isidentical(Z2, Zonotope(N[4, 5, 6], N[3 7; 3 7; 3 7]))

    # distance (between point and set)
    @test_broken distance(Z, N[0, 0]) isa LazySet{N}  # TODO implement `distance` for polytopes
    # @test_throws DimensionMismatch distance(Z, N[0])
    # @test_throws ArgumentError distance(Z, N[0]; p=N(1 // 2))
    # for (x, v) in ((N[1], N(0)), (N[4], N(2)))
    #     for res in (distance(Z, x), distance(x, Z))
    #         @test res == v
    #         if N <: AbstractFloat
    #             @test res isa N
    #         end
    #     end
    # end

    # exponential_map (part 1)
    @test_throws DimensionMismatch exponential_map(ones(N, 1, 1), Z)
    @test_throws DimensionMismatch exponential_map(ones(N, 2, 1), Z)

    # in
    @test_throws DimensionMismatch N[0] ∈ Z
    @test (@inferred N[1, 2] ∈ Z) && (@inferred N[5, 8] ∈ Z) && @inferred N[3, 3] ∉ Z
    @test (@inferred N[1, 2] ∈ Z0) && @inferred N[2, 3] ∉ Z0

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0], Z)
    @test_throws ArgumentError is_interior_point(N[0, 0], Z; ε=N(0))
    @test_throws ArgumentError is_interior_point(N[0, 0], Z; p=N(1 // 2))
    if N <: AbstractFloat
        @test is_interior_point(N[1, 2], Z)
        @test_broken @inferred is_interior_point(N[5, 8], Z)  # TODO make this type-stable (witness)
        @test !is_interior_point(N[5, 8], Z)
        @test !is_interior_point(N[3, 3], Z)
    else
        @test_throws ArgumentError is_interior_point(N[0, 0], Z)
        @test is_interior_point(N[1, 2], Z; ε=1 // 100)
        @test_broken !(@inferred is_interior_point(N[5, 8], Z; ε=1 // 100))  # TODO make this type-stable (witness)
        @test !is_interior_point(N[5, 8], Z; ε=1 // 100)
        @test !is_interior_point(N[3, 3], Z; ε=1 // 100)
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([0.0, 0.0], Z)
    end

    # linear_map
    @test_throws DimensionMismatch linear_map(ones(N, 1, 1), Z)
    Z2 = @inferred linear_map(N[1 0; 0 2], Z)
    @test isidentical(Z2, Zonotope(N[1, 4], N[1 3; 4 8]))
    Z2 = @inferred linear_map(zeros(N, 1, 2), Z)  # zero map
    @test isidentical(Z2, Zonotope(N[0], zeros(N, 1, 0)))
    # 1D simplifies to 1 generator (except for the zero map)
    Z2 = @inferred linear_map(N[-1 1;], Z)
    @test isidentical(Z2, Zonotope(N[1], hcat(N[2])))
    # linear_map!
    Z2 = Zonotope(similar(Z.center), similar(Z.generators))
    @inferred linear_map!(Z2, N[1 0; 0 2], Z)
    @test isidentical(Z2, Zonotope(N[1, 4], N[1 3; 4 8]))

    # linear_map_inverse
    @test_throws AssertionError LazySets.linear_map_inverse(ones(N, 1, 1), Z)
    # invertible map
    Z2 = LazySets.linear_map_inverse(N[1 0; 0 1], Z)
    @test isidentical(Z2, Z)
    # noninvertible map
    M = ones(N, 2, 1)
    X = LazySets.linear_map_inverse(M, Z)
    A, b = tosimplehrep(Z)
    Y = HPolytope(A * M, b)
    @test X isa LazySet{N} && isequivalent(X, Y)

    # permute
    @test_throws DimensionMismatch permute(Z, [1])
    @test_throws DimensionMismatch permute(Z, [1, -1])
    @test_throws DimensionMismatch permute(Z, [1, 3])
    @test_throws ArgumentError permute(Z, [1, 1])
    Z2 = @inferred permute(Z, [1, 2])
    @test isidentical(Z2, Z)
    Z2 = @inferred permute(Z, [2, 1])
    @test isidentical(Z2, Zonotope(N[2, 1], N[2 4; 1 3]))

    # project
    @test_throws DimensionMismatch project(Z, [1, 2, 3])
    @test_throws DimensionMismatch project(Z, [1, -1])
    @test_throws DimensionMismatch project(Z, [1, 3])
    @test_throws ArgumentError project(Z, [1, 1])
    Z2 = @inferred project(Z, [1])
    @test isidentical(Z2, Zonotope(N[1], hcat(N[4])))
    # 1D projection with zero generators (#2147)
    Z2 = @inferred project(convert(Zonotope, BallInf(N[0, 0], N(1))), [1])
    @test isidentical(Z2, Zonotope(N[0], hcat(N[1])))
    # 2D projection
    Z2a = Zonotope(zeros(N, 3), N[1 0 3; 4 0 6; 0 0 0])
    Z2 = @inferred project(Z2a, [1, 2])
    # removing zero generators in projection can be deactivated
    @test isidentical(Z2, Zonotope(zeros(N, 2), N[1 3; 4 6]))
    Z2 = @inferred project(Z2a, [1, 2]; remove_zero_generators=false)
    @test isidentical(Z2, Zonotope(zeros(N, 2), N[1 0 3; 4 0 6]))

    # sample
    res = @inferred sample(Z)
    @test res isa Vector{N} && res ∈ Z
    res = @inferred sample(Z, 2)
    @test res isa Vector{Vector{N}} && length(res) == 2 && all(x ∈ Z for x in res)

    # scale
    Z2 = @inferred scale(N(2), Z)
    @test isidentical(Z2, Zonotope(N[2, 4], N[2 6; 4 8]))
    # degenerate case
    Z2 = @inferred scale(N(0), Z)
    @test isidentical(Z2, Zonotope(N[0, 0], zeros(N, 2, 0)))
    # scale!
    Z2 = copy(Z)
    @inferred scale!(N(2), Z2)
    @test isidentical(Z2, Zonotope(N[2, 4], N[2 6; 4 8]))
    # degenerate case
    Z2 = copy(Z)
    @inferred scale!(N(0), Z2)
    @test isidentical(Z2, Zonotope(N[0, 0], zeros(N, 2, 2)))

    # split
    # split into two zonotopes
    @test_throws AssertionError split(Z, 0)
    Z2a, Z2b = @inferred split(Z, 2)
    @test isidentical(Z2a, Zonotope(N[-1 // 2, 0], N[1 3//2; 2 2]))
    @test isidentical(Z2b, Zonotope(N[5 // 2, 4], N[1 3//2; 2 2]))
    # split into multiple zonotopes
    Z2a = Zonotope(N[0, 0], N[1 1; -1 1])
    Z2a, Z2b, Z2c, Z2d = @inferred split(Z2a, [1, 2], [1, 1])
    G = N[1//2 1//2; -1//2 1//2]
    ispermutation([Z2a, Z2b, Z2c, Z2d],
                  [Zonotope(N[-1, 0], G), Zonotope(N[0, 1], G),
                   Zonotope(N[0, -1], G), Zonotope(N[1, 0], G)])
    @static if isdefined(@__MODULE__, :StaticArrays)
        using StaticArrays: SVector, SMatrix

        Z2 = Zonotope(SVector{2}(N[0, 0]), SMatrix{2,2}(N[1 1; -1 1]))
        Z2a, Z2b = @inferred split(Z2, 1)
        @test Z2a ⊆ Z2 && Z2b ⊆ Z2
    end

    # support_function
    @test_throws DimensionMismatch ρ(N[1], Z)
    res = @inferred ρ(N[1, 1], Z)
    @test res isa N && res == N(13)
    res = @inferred ρ(N[-1, -1], Z)
    @test res isa N && res == N(7)

    # support_vector
    @test_throws DimensionMismatch σ(N[1], Z)
    res = @inferred σ(N[1, 1], Z)
    @test res isa Vector{N} && res == N[5, 8]
    res = @inferred σ(N[-1, -1], Z)
    @test res isa Vector{N} && res == N[-3, -4]

    # translate
    @test_throws DimensionMismatch translate(Z, N[1])
    Z2 = @inferred translate(Z, N[1, 2])
    @test isidentical(Z2, Zonotope(N[2, 4], N[1 3; 2 4]))
    # translate!
    @test_throws DimensionMismatch translate!(Z, N[1])
    Z2 = copy(Z)
    @inferred translate!(Z2, N[1, 2])
    @test isidentical(Z2, Zonotope(N[2, 4], N[1 3; 2 4]))
    # test immutable usage
    @static if isdefined(@__MODULE__, :StaticArrays)
        using StaticArrays: SVector, SMatrix

        Z2 = Zonotope(SVector{2}(c), SMatrix{2,2}(N[1 3; 2 4]))
        Z2b = @inferred translate(Z2, N[1, 2])
        @test isidentical(Z2b, Zonotope(SVector{2}(N[2, 4]), SMatrix{2,2}(N[1 3; 2 4])))
        @test_throws ErrorException translate!(Z2, N[1, 2])
    end

    # cartesian_product
    Z2 = Zonotope(N[3], hcat(N[1]))
    Z2b = @inferred cartesian_product(Z, Z2)
    @test isidentical(Z2b, Zonotope(N[1, 2, 3], N[1 3 0; 2 4 0; 0 0 1]))
    Z2b = @inferred cartesian_product(Z2, Z)
    @test isidentical(Z2b, Zonotope(N[3, 1, 2], N[1 0 0; 0 1 3; 0 2 4]))

    # convex_hull (binary)
    @test_throws DimensionMismatch convex_hull(Z, Z3)
    X = convex_hull(Z, Z)
    @test X isa LazySet{N} && isequivalent(X, Z)
    Z2 = Zonotope(N[4, -2], N[1 0; 0 2])
    Y = VPolygon([N[-3, -4], N[5, -4], N[5, 8], N[-1, 0]])
    for X in (convex_hull(Z, Z2), convex_hull(Z2, Z))
        @test X isa LazySet{N} && isequivalent(X, Y)
    end

    # difference
    @test_broken difference(Z, Z3) isa LazySet{N}  # TODO implement `difference` for polytopes (add more tests later)

    # distance (between two sets)
    @test_broken distance(Z, Z3) isa LazySet{N}  # TODO implement `distance` for polytopes (add more tests later)

    # exact_sum
    @test_throws DimensionMismatch exact_sum(Z, Z3)
    Z2 = Zonotope(N[4, -2], N[1 0; 0 2])
    for Z2b in ((@inferred exact_sum(Z, Z2)), @inferred exact_sum(Z2, Z))
        @test Z2b isa Zonotope{N} && isequivalent(Z2b, Zonotope(N[5, 0], N[1 3 1 0; 2 4 0 2]))
    end

    # intersection
    @test_throws DimensionMismatch intersection(Z, Z3)
    # disjoint
    X = intersection(Z, Zonotope(N[3, 0], N[2 0; 0 1]))
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    # overlapping
    X = intersection(Z, Zonotope(N[1, 2], N[2 0; 0 1]))
    Y = VPolygon([N[5 // 4, 3], N[-1 // 4, 1], N[3 // 4, 1], N[9 // 4, 3]])
    if N != Float32
        @test X isa LazySet{N} && isequivalent(X, Y)
    end

    # isapprox
    @test @inferred Z ≈ Z
    res = @inferred Z ≈ translate(Z, N[1 // 100000000, 0])
    if N <: AbstractFloat
        @test res  # below default tolerance for AbstractFloat
    else
        @test !res  # zero default tolerance for Rational
    end
    @test !(@inferred Z ≈ translate(Z, N[1 // 1000, 0]))  # above default tolerance for all types
    @test !(@inferred Z ≈ Z3) && !(@inferred Z3 ≈ Z) && !(@inferred Z ≈ P) && !(@inferred P ≈ Z)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(Z, Z3)
    # disjoint
    Z2 = Zonotope(N[2, -2], N[1 0; 0 1])
    @test_broken @inferred isdisjoint(Z, Z2)  # TODO make this type-stable (witness)
    @test isdisjoint(Z, Z2) && isdisjoint(Z2, Z)
    for (Z2a, Z2b) in ((Z, Z2), (Z2, Z))
        @test_broken @inferred isdisjoint(Z2a, Z2b, true)  # TODO make this type-stable (witness)
        res, w = isdisjoint(Z2a, Z2b, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    # overlapping
    Z2 = Zonotope(N[2, 2], N[1 0 1; 0 1 1])
    @test !isdisjoint(Z, Z) && !isdisjoint(Z, Z2) && !isdisjoint(Z2, Z)
    for (Z2a, Z2b) in ((Z, Z), (Z, Z2), (Z2, Z))
        res, w = isdisjoint(Z2a, Z2b, true)
        @test !res && w isa Vector{N} && w ∈ Z2a && w ∈ Z2b
    end
    # no generators
    Z2 = Zonotope(N[2, 3], zeros(N, 2, 0))
    @test isdisjoint(Z0, Z2) && isdisjoint(Z2, Z0)
    for (Z2a, Z2b) in ((Z0, Z2), (Z2, Z0))
        res, w = isdisjoint(Z2a, Z2b, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    @test !isdisjoint(Z0, Z0)
    res, w = isdisjoint(Z0, Z0, true)
    @test !res && w isa Vector{N} && w ∈ Z0
    # tolerance
    if N == Float64
        Z2 = Zonotope(N[6 + 1e-9, 8], N[1 0; 0 1])
        @test !isdisjoint(Z, Z2)
        r = _rtol(N)
        @assert r > N(1e-10) "default tolerance changed; adapt test"
        set_rtol(N, N(1e-10))
        @test_broken isdisjoint(Z, Z2)  # cannot adapt tolerance of LP solver
        # restore tolerance
        set_rtol(N, r)
    end

    # isequal
    @test @inferred Z == Z
    @test (@inferred Z != Z3) && (@inferred Z3 != Z) && (@inferred Z != P) && @inferred P != Z

    # isequivalent
    @test_throws DimensionMismatch isequivalent(Z, Z3)
    @test_broken @inferred isequivalent(Z, Z)  # TODO make this type-stable (witness)
    @test isequivalent(Z, Z)
    @test !isequivalent(Z, Zonotope(N[1, 2], N[1 0; 0 1]))
    @test isequivalent(Z, P) && isequivalent(P, Z)

    # isstrictsubset
    @test_throws DimensionMismatch Z ⊂ Z3
    Z2 = Zonotope(N[1, 2], N[1 0; 0 1])
    @test_broken !(@inferred Z2 ⊂ Z)  # TODO make this type-stable (witness)
    @test !(Z2 ⊂ Z)
    @test_broken @inferred ⊂(Z2, Z, true)  # TODO make this type-stable (witness)
    res, w = ⊂(Z2, Z, true)
    @test !res && w isa Vector{N} && w ∈ Z2 && w ∉ Z
    for Z2 in (Z, P)
        @test !(Z2 ⊂ Z)
        res, w = ⊂(Z2, Z, true)
        @test !res && w isa Vector{N} && isempty(w)
    end
    Z2 = Zonotope(N[1, 2], N[4 0; 0 6])
    @test Z ⊂ Z2
    res, w = ⊂(Z, Z2, true)
    @test res && w isa Vector{N} && w ∉ Z && w ∈ Z2

    # issubset
    @test_throws DimensionMismatch Z ⊆ Z3
    for X in (Z, P)
        @test_broken @inferred Z ⊆ X  # TODO make this type-stable (witness)
        @test Z ⊆ X
        @test_broken @inferred ⊆(Z, X, true)  # TODO make this type-stable (witness)
        res, w = ⊆(Z, X, true)
        @test res && w isa Vector{N} && w == N[]
    end
    Z2 = Zonotope(N[3, 0], N[1 0; 0 1])
    @test Z ⊈ Z2
    res, w = ⊆(Z, Z2, true)
    @test !res && w isa Vector{N} && w ∈ Z && w ∉ Z2

    # linear_combination
    @test_throws DimensionMismatch linear_combination(Z, Z3)
    @test_broken linear_combination(Z, Xnc) isa LazySet{N}  # TODO implement `linear_combination` for non-convex sets
    @test_broken linear_combination(Xnc, Z) isa LazySet{N}
    for X in (linear_combination(Z, Z), linear_combination(Z, P), linear_combination(P, Z))
        @test X isa LazySet{N} && isequivalent(X, Z)
    end
    Z2 = Zonotope(N[4, -2], N[1 0; 0 2])
    Y = VPolygon([N[-3, -4], N[5, -4], N[5, 8], N[-1, 0]])
    for X in (linear_combination(Z, Z2), linear_combination(Z2, Z))
        @test X isa LazySet{N} && isequivalent(X, Y)
    end

    # minkowski_difference (part 1)
    # equivalent sets
    if VERSION >= v"1.11"
        for X in ((@inferred minkowski_difference(Z, P)), @inferred minkowski_difference(P, Z))
            vlist = vertices_list(X)
            @test length(vlist) == 1 && all(isapproxzero, vlist[1])
        end
    else
        for X in (minkowski_difference(Z, P), minkowski_difference(P, Z))
            vlist = vertices_list(X)
            @test length(vlist) == 1 && all(isapproxzero, vlist[1])
        end
    end

    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(Z, Z3)
    Z2 = @inferred minkowski_sum(Z, Z)
    @test isidentical(Z2, Zonotope(N[2, 4], N[1 3 1 3; 2 4 2 4]))
end

for N in @tN([Float64, Float32])
    Z = Zonotope(N[1, 2], N[1 3; 2 4])
    Z3 = Zonotope(N[1, 1, 1], N[1 0 0; 0 2 0; 0 0 3])

    # rand
    Z2 = rand(Zonotope; N=N)
    @test Z2 isa Zonotope{N} && dim(Z2) == 2
    Z2 = rand(Zonotope; N=N, dim=3)
    @test Z2 isa Zonotope{N} && dim(Z2) == 3
    Z2 = rand(Zonotope; N=N, num_generators=5)
    @test Z2 isa Zonotope{N} && dim(Z2) == 2 && ngens(Z2) == 5

    # rationalize
    @test_broken @inferred rationalize(Z)  # TODO make this type-stable (fix in ReachabilityBase)
    Z2 = rationalize(Z)
    T = Rational{Int64}
    @test Z2 isa Zonotope{T} && isidentical(Z2, Zonotope(T[1, 2], T[1 3; 2 4]))

    # exponential_map (part 2)
    Z2 = @inferred exponential_map(N[1 0; 0 2], Z)
    @test isidentical(Z2, Zonotope(N[exp(1), 2 * exp(2)], N[exp(1) 3*exp(1); 2*exp(2) 4*exp(2)]))

    # minkowski_difference (part 2)
    @test_throws DimensionMismatch minkowski_difference(Z, Z3)
    # empty difference
    Z2 = Zonotope(N[1, 2], N[4 0; 0 6])
    X = minkowski_difference(Z, Z2)
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    # nonempty difference
    Z2 = minkowski_difference(Z, Z)
    @test isidentical(Z2, Zonotope(N[0, 0], zeros(N, 2, 0)))
end

for N in @tN([Float64, Rational{Int}])
    Z = Zonotope(N[1, 2], N[1 3; 2 4])
    Z3 = Zonotope(N[1, 1, 1], N[1 0 0; 0 2 0; 0 0 3])  # 3D zonotope, a hyperrectangle

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        X = polyhedron(Z)
        @test X isa Polyhedra.DefaultPolyhedron{N}
        @test ispermutation(constraints_list(convert(HPolytope, X)), constraints_list(Z))
    end

    # triangulate
    @static if isdefined(@__MODULE__, :MiniQhull)
        @static if VERSION >= v"1.11"
            X = @inferred triangulate(Z)
        else
            X = triangulate(Z)
        end
        @test X isa UnionSetArray{N,<:VPolytope{N}} && length(X) == 2
        # result is not unique; reylying on implementation here
        for Xi in X
            @test ispermutation(vertices_list(Xi), [N[3, 4], N[-1, 0], N[-3, -4]]) ||
                  ispermutation(vertices_list(Xi), [N[3, 4], N[-1, 0], N[5, 8]])
        end
    end

    # volume (part 2)
    @static if isdefined(@__MODULE__, :Polyhedra)
        if N <: AbstractFloat && VERSION >= v"1.11"
            res = @inferred volume(Z3)
        else
            res = volume(Z3)
        end
        @test res isa N && res == N(48)
    end
end

for N in [Float64]
    Z3 = Zonotope(N[1, 1, 1], N[1 0 0; 0 2 0; 0 0 3])  # 3D zonotope, a hyperrectangle

    # chebyshev_center_radius
    @static if isdefined(@__MODULE__, :Polyhedra)
        @static if VERSION >= v"1.12"
            c, r = @inferred chebyshev_center_radius(Z3)
        else
            c, r = chebyshev_center_radius(Z3)
        end
        @test c isa Vector{N} && c == Z3.center
        @test r isa N && r == N(1)
    end
end
