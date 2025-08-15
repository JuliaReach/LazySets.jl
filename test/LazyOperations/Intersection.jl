using LazySets, Test
using LazySets.ReachabilityBase.Arrays: SingleEntryVector, ispermutation
@static if VERSION >= v"1.9"
    vGLPK = pkgversion(LazySets.GLPK)
else
    import PkgVersion
    vGLPK = PkgVersion.Version(LazySets.GLPK)
end
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    B = BallInf(ones(N, 2), N(3))
    H = Hyperrectangle(ones(N, 2), ones(N, 2))
    E = EmptySet{N}(2)

    # intersection of two sets
    I = Intersection(B, H)

    # emptiness of intersection
    @test !isempty_known(I)
    @test !isempty(I)
    @test isempty_known(I)
    @test !isempty(I)

    # convenience constructors
    @test B ∩ H == I
    cap = IntersectionArray([B, H, B])
    @test ∩(B, H, B) == ∩([B, H, B]) == cap
    @test ∩(B) == B

    # Universe is neutral
    U = Universe{N}(2)
    @test B ∩ U == U ∩ B == B
    @test U ∩ U == U
    ia = IntersectionArray([B, N(2) * B, N(3) * B])
    @test ia ∩ U == U ∩ ia == ia

    # array interface
    @test array(I) == [B, H] && array(cap) == [B, H, B]
    @test I[1] == cap[1] == B
    @test I[1:2] == cap[1:2] == [B, H]
    @test I[end] == H && cap[end] == B
    @test length(I) == 2 && length(cap) == 3
    v = Vector{LazySet{N}}()
    @test array(IntersectionArray(v)) ≡ v

    # swap
    I2 = swap(I)
    @test I.X == I2.Y && I.Y == I2.X

    # flatten
    B2 = Ball1(N[0, 0], N(1))
    for M3 in (Intersection(Intersection(B, IntersectionArray([H])), B2),
               IntersectionArray([Intersection(B, IntersectionArray([H])), B2]))
        M3f = flatten(M3)
        @test M3f isa IntersectionArray && array(M3f) == [B, H, B2]
    end

    # dim
    @test dim(I) == 2

    # support vector
    @test σ(ones(N, 2), I) == N[2, 2]

    # membership
    @test ones(N, 2) ∈ I && N[5, 5] ∉ I

    # boundedness (more tests in Float64-only section)
    @test isbounded(I) && isboundedtype(typeof(I))
    I2 = Singleton(N[1]) ∩ HalfSpace(N[1], N(1))
    @test isbounded(I2) && isboundedtype(typeof(I2))
    I2 = Hyperplane(ones(N, 2), N(1)) ∩ HalfSpace(ones(N, 2), N(1))
    # the next test fails with the "EXACT" solver GLPK older than v0.15.3
    # see https://github.com/jump-dev/GLPK.jl/issues/207
    if N != Rational{Int} || vGLPK >= v"0.15.3"
        @test !isbounded(I2) && !isboundedtype(typeof(I2))
    end

    # ispolyhedral
    @test ispolyhedral(I)
    if N isa AbstractFloat
        I2 = Intersection(B, Ball2(N[0, 0], N(1)))
        @test !ispolyhedral(I2)
    end

    # concretize
    @test LazySets.concrete_function(Intersection) == intersection
    @test concretize(I) == intersection(B, H)

    # volume
    @test volume(I) == N(4)

    # conversion of intersection of polyhedral types to polyhedral types
    H1 = HalfSpace(N[1, 1], N(1))
    H2 = HalfSpace(N[-1, -1], N(1))
    P12 = convert(HPolyhedron, H1 ∩ H2)
    @test P12 == HPolyhedron([H1, H2])

    # simplification of lazy intersection of half-spaces to polyhedron
    @test H1 ∩ H2 == HPolyhedron([H1, H2])
    H3 = HalfSpace(N[1, 0], N(3))
    @test H1 ∩ H2 ∩ H3 == HPolyhedron([H1, H2, H3])
    # make lazy intersection of half-spaces by specifying the cache
    Ih12 = Intersection(H1, H2; cache=LazySets.IntersectionCache())
    @test Ih12.X == H1 && Ih12.Y == H2

    # vertices
    I2 = Intersection(BallInf(zeros(N, 2), N(2)), BallInf(3 * ones(N, 2), N(2)))
    @test ispermutation(vertices_list(I2),
                        vertices_list(BallInf(fill(N(3 // 2), 2), N(1 // 2))))

    # =================
    # IntersectionArray
    # =================

    # relation to base type (internal helper functions)
    @test LazySets.array_constructor(Intersection) == IntersectionArray
    @test LazySets.binary_constructor(IntersectionArray) == Intersection
    @test !LazySets.is_array_constructor(Intersection)
    @test LazySets.is_array_constructor(IntersectionArray)

    # intersection of an array of sets
    IArr = IntersectionArray([B, H])

    # getindex & length
    @test IArr[1] == B && IArr[2] == H
    @test length(IArr) == 2

    # dim
    @test dim(IArr) == 2

    # support vector
    @test σ(ones(N, 2), IArr) == N[2, 2]

    # boundedness
    @test isbounded(IArr) && isboundedtype(typeof(IArr))
    @test isbounded(IntersectionArray([Singleton(N[1]), HalfSpace(N[1], N(1))]))
    @test isbounded(IntersectionArray([HalfSpace(N[1], N(1)), HalfSpace(N[-1], N(1))]))
    IArr2 = IntersectionArray([HalfSpace(N[1], N(1))])
    @test !isboundedtype(typeof(IArr2))

    # ispolyhedral
    @test ispolyhedral(IArr)
    if N isa AbstractFloat
        IArr2 = IntersectionArray([B, Ball2(N[0, 0], N(1))])
        @test !ispolyhedral(IArr2)
    end

    # isempty
    @test !isempty(IArr)

    # membership
    @test ones(N, 2) ∈ IArr && N[5, 5] ∉ IArr

    # constructor with size hint and type
    IntersectionArray(10, N)

    # concretize
    @test concretize(IArr) == intersection(B, H)

    # conversion of intersection of polyhedral types to polyhedral types
    P12 = convert(HPolyhedron, IntersectionArray([H1, H2]))
    @test P12 == HPolyhedron([H1, H2])

    # ================
    # common functions
    # ================

    # absorbing element
    @test absorbing(Intersection) == absorbing(IntersectionArray) == EmptySet
    @test I ∩ E == E ∩ I == IArr ∩ E == E ∩ IArr == E ∩ E == E

    # intersection between a hyperrectangular set and an axis-aligned halfspace
    SV = SingleEntryVector(1, 2, N(1))
    H = HalfSpace(SV, N(2))
    B1 = Hyperrectangle(N[1, 1], N[1, 1])
    @test intersection(H, B1) == B1
    B2 = Hyperrectangle(N[2, 1], N[1, 1])
    @test intersection(H, B2) == Hyperrectangle(N[1.5, 1], N[0.5, 1])
    B3 = Hyperrectangle(N[10, 10], N[1, 1])
    @test intersection(H, B3) == EmptySet{N}(2)
end

for N in [Float64]
    # constraints_list for polytopic intersection
    B = BallInf(ones(N, 2), N(3))
    H = Hyperrectangle(ones(N, 2), ones(N, 2))
    I = B ∩ H
    IArr = IntersectionArray([B, H])
    clist1 = constraints_list(I)
    clist2 = constraints_list(IArr)
    @test ispermutation(clist1, clist2) &&
          ispermutation(clist1,
                        [HalfSpace(N[1, 0], N(2)),
                         HalfSpace(N[0, 1], N(2)),
                         HalfSpace(N[-1, 0], N(0)),
                         HalfSpace(N[0, -1], N(0))])

    # constraints_list with single HalfSpace
    H2 = HalfSpace(N[1], N(0))
    IArr = IntersectionArray([H2])
    @test constraints_list(IArr) == [H2]

    # HalfSpace vs. Ball1 intersection
    X = Ball1(zeros(2), N(1))
    d = normalize(N[1, 0])

    # flat intersection at x = 1
    H = HalfSpace(N[-1, 0], N(-1)) # x >= 1

    # default algorithm
    @static if isdefined(@__MODULE__, :Optim)
        @test ρ(d, X ∩ H) == ρ(d, H ∩ X) == N(1)
    end

    # intersection at x = 0
    H = HalfSpace(N[1, 0], N(0)) # x <= 0

    # default algorithm
    @static if isdefined(@__MODULE__, :Optim)
        @test ρ(d, X ∩ H) < 1e-6 && ρ(d, H ∩ X) < 1e-6
    end

    # specify line-search algorithm
    @static if isdefined(@__MODULE__, :Optim)
        @test ρ(d, X ∩ H; algorithm="line_search") < 1e-6 &&
              ρ(d, H ∩ X; algorithm="line_search") < 1e-6

        # test approximation errors (see #1144)
        B = Hyperrectangle(N[0.009231278330571413, 0], N[0.009231221669425305, 1])
        H = HalfSpace(N[1, 0], N(0))
        @test ρ([1.0, 0], B ∩ H) >= 0
    end

    # boundedness
    @test isbounded(HalfSpace(N[1], N(1)) ∩ HalfSpace(N[-1], N(1)))
    @test !isbounded(HalfSpace(ones(N, 2), N(1)) ∩ HalfSpace(-ones(N, 2), N(1)))
    @test !isbounded(IntersectionArray([HalfSpace(N[1], N(1)), HalfSpace(N[1], N(-1))]))
    @test !isbounded(IntersectionArray([HalfSpace(ones(N, 2), N(1)), HalfSpace(ones(N, 2), N(-1))]))

    # HalfSpace vs. Ball2 intersection
    @static if isdefined(@__MODULE__, :Optim)
        B2 = Ball2(zeros(2), N(1))
        @test ρ(d, B2 ∩ H) < 1e-6 && ρ(d, H ∩ B2) < 1e-6

        # Ball1 vs. Hyperplane intersection
        H = Hyperplane(N[1, 0], N(0.5)) # x = 0.5
        @test isapprox(ρ(d, X ∩ H; algorithm="line_search"), N(0.5), atol=1e-6)
        # For the projection algorithm, if the linear map is taken lazily we can use Ball1
        @test isapprox(ρ(d, X ∩ H; algorithm="projection", lazy_linear_map=true), N(0.5), atol=1e-6)
        # But the default is to take the linear map concretely; in this case, we *may*
        # need Polyhedra (in the general case), for the concrete linear map. As a valid workaround
        # if we don't want to load Polyhedra here is to convert the given set to a polygon in V-representation
        @test isapprox(ρ(d, convert(VPolygon, X) ∩ H; algorithm="projection",
                         lazy_linear_map=false), N(0.5), atol=1e-6)
    end

    # =====================
    # concrete operations
    # =====================
    cap = HPolytope([HalfSpace(N[1], N(1))]) ∩ HPolytope([HalfSpace(N[-1], N(1))])  # x <= 1 && x >= -1
    p = linear_map(reshape([N(1 / 2)], 1, 1), cap)
    @test N[-0.5] ∈ p && N[0.5] ∈ p && (N[1.0] ∉ p || N[1.0] ∉ p)
end
