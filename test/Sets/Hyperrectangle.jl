using LazySets, Test, SparseArrays, LinearAlgebra
using LazySets.ReachabilityBase.Arrays: ispermutation, SingleEntryVector
using LazySets.ReachabilityBase.Comparison: set_rtol, _rtol
IA = LazySets.IA
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

function isidentical(::Hyperrectangle, ::Hyperrectangle)
    return false
end

function isidentical(H1::Hyperrectangle{N}, H2::Hyperrectangle{N}) where {N}
    return H1.center == H2.center && H1.radius == H2.radius
end

for N in @tN([Float64, Float32, Rational{Int}])
    # auxiliary sets
    P = VPolygon([N[2, 1], N[0, 1], N[0, -3], N[2, -3]])  # equivalent set
    Xnc = UnionSet(P, BallInf(N[3, 0], N(1)))  # nonconvex set

    # constructor
    c = N[1, -1]
    r = N[1, 2]
    @test_throws AssertionError Hyperrectangle(c, N[1])
    @test_throws AssertionError Hyperrectangle(c, N[1, -1])
    @test_throws AssertionError Hyperrectangle(c, N[1, Inf])
    H = @inferred Hyperrectangle(c, r)
    H1 = @inferred Hyperrectangle(N[0], N[1])  # 1D hyperrectangle
    H3 = @inferred Hyperrectangle(N[0, 0, 0], N[1, 2, 3])  # 3D hyperrectangle
    H0 = @inferred Hyperrectangle(N[0, 0], N[0, 0])  # flat hyperrectangle
    # constructor from low/high vectors
    @test_throws AssertionError Hyperrectangle(low=N[1], high=N[0])
    @test (@inferred Hyperrectangle(; low=N[1], high=N[0], check_bounds=false)) isa
          Hyperrectangle{N}
    H2 = @inferred Hyperrectangle(; low=c - r, high=c + r)
    @test H2 == H
    # constructor with mixed vectors
    H2 = @inferred Hyperrectangle(sparsevec([1], N[0], 1), N[1])
    @test H2 == H1
    H2 = @inferred Hyperrectangle(N[0], sparsevec([1], N[1], 1))
    @test H2 == H1
    # unicode constructor
    @test isidentical(□(c, r), H)

    # convert
    # from AbstractHyperrectangle
    H2 = @inferred convert(Hyperrectangle, BallInf(N[2, 3], N(1)))
    @test isidentical(H2, Hyperrectangle(N[2, 3], N[1, 1]))
    # between IntervalBox
    @static if isdefined(@__MODULE__, :IntervalBoxes)
        import IntervalBoxes as IB

        @test_broken @inferred convert(IB.IntervalBox, H)  # TODO make this type-stable
        X = convert(IB.IntervalBox, H)
        @test X isa IB.IntervalBox{2,N} && Interval(X[1]) == Interval(N(0), N(2)) &&
              Interval(X[2]) == Interval(N(-3), N(1))
        @test isidentical(convert(Hyperrectangle, X), H)
    end
    # from IA.Interval
    H2 = @inferred convert(Hyperrectangle, IA.interval(N(-1), N(1)))
    @test isidentical(H2, H1)
    # from CartesianProduct of AbstractHyperrectangles
    X = Interval(N(0), N(2)) × Interval(N(-3), N(1))
    @test isidentical((@inferred convert(Hyperrectangle, X)), H)
    # from CartesianProductArray of Intervals
    X = CartesianProductArray([X[1], X[2]])
    @test isidentical((@inferred convert(Hyperrectangle, X)), H)
    # from CartesianProductArray of AbstractHyperrectangles
    X = CartesianProductArray(convert.(Hyperrectangle, array(X)))
    @test isidentical((@inferred convert(Hyperrectangle, X)), H)
    # from Rectification of AbstractHyperrectangle
    H2 = @inferred convert(Hyperrectangle, Rectification(H))
    @test isidentical(H2, Hyperrectangle(N[1, 1 // 2], N[1, 1 // 2]))
    # from Zonotope
    X = convert(Zonotope, H)
    @test isidentical((@inferred convert(Hyperrectangle, X)), H)
    X = Zonotope(c, N[2 1; 1 2])
    @test_throws AssertionError convert(Hyperrectangle, X)
    X = Zonotope(c, N[1 0 0; 0 0 2])  # zero generator
    @test isidentical((@inferred convert(Hyperrectangle, X)), H)

    # an_element
    x = @inferred an_element(H)
    @test x isa Vector{N} && x ∈ H

    # area
    @test_throws DimensionMismatch area(H1)
    res = @inferred area(H)
    @test res isa N && res == N(8)
    res = @inferred area(H0)
    @test res isa N && res == N(0)
    res = @inferred area(H3)
    @test res isa N && res == N(88)

    # center
    c2 = @inferred center(H)
    @test c2 isa AbstractVector{N} && c2 == c
    @test_throws DimensionMismatch center(H, 3)
    v = @inferred center(H, 1)
    @test v isa N && v == N(1)

    # chebyshev_center_radius
    c2, r2 = @inferred chebyshev_center_radius(H)
    @test c2 isa Vector{N} && c2 == c && r2 isa N && r2 == N(1)

    # complement
    X = @inferred complement(H)
    @test X isa UnionSetArray{N} && dim(X) == 2
    clist = [HalfSpace(N[-1, 0], N(-2)), HalfSpace(N[0, -1], N(-1)),
             HalfSpace(N[1, 0], N(0)), HalfSpace(N[0, 1], N(-3))]
    @test ispermutation(array(X), clist)

    # concretize
    H2 = @inferred concretize(H)
    @test isidentical(H2, H)

    # constrained_dimensions
    @test (@inferred constrained_dimensions(H)) == 1:2

    # constraints_list
    clist = [HalfSpace(N[1, 0], N(2)), HalfSpace(N[0, 1], N(1)),
             HalfSpace(N[-1, 0], N(0)), HalfSpace(N[0, -1], N(3))]
    @test ispermutation((@inferred constraints_list(H)), clist)

    # constraints
    @test ispermutation(collect(@inferred constraints(H)), clist)

    # convex_hull (unary)
    H2 = @inferred convex_hull(H)
    @test isidentical(H, H2)

    # copy
    H2 = @inferred copy(H)
    @test isidentical(H, H2)

    # diameter
    @test_throws ArgumentError diameter(H, N(1 // 2))
    for res in ((@inferred diameter(H)), @inferred diameter(H, Inf))
        if N <: AbstractFloat
            @test res isa N && res == N(4)
        else
            @test res == 4.0
        end
    end
    res = @inferred diameter(H, 2)
    if N <: AbstractFloat
        @test res isa N
    end
    @test res == 2 * norm(N[1, 2], 2)

    # dim
    @test (@inferred dim(H)) == 2
    @test (@inferred dim(H3)) == 3

    # eltype
    @test (@inferred eltype(H)) == N
    @test (@inferred eltype(typeof(H))) == N

    # extrema
    res = @inferred extrema(H)
    @test res isa Tuple{Vector{N},Vector{N}} && res[1] == N[0, -3] && res[2] == N[2, 1]
    @test_throws DimensionMismatch extrema(H, 3)
    res = @inferred extrema(H, 1)
    @test res isa Tuple{N,N} && res[1] == N(0) && res[2] == N(2)

    # generators
    @test ispermutation(collect((@inferred generators(H))), [N[1, 0], N[0, 2]])
    # degenerate case
    gens = collect(@inferred generators(H0))
    @test gens isa AbstractVector{<:AbstractVector{N}} && isempty(gens)

    # genmat
    @test (@inferred genmat(H)) == N[1 0; 0 2]
    # degenerate case
    gens = @inferred genmat(H0)
    @test gens isa Matrix{N} && gens == Matrix{N}(undef, 2, 0)
    # genmat with sparse or static arrays
    @static if isdefined(@__MODULE__, :StaticArrays)
        using StaticArrays: SA, SMatrix

        G1 = @inferred genmat(Hyperrectangle(sparsevec(N[3, 2]), sparsevec(N[2, 1])))
        @test_broken @inferred genmat(Hyperrectangle(SA[N(3), N(2)], SA[N(2), N(1)]))  # TODO make this type-stable
        G2 = genmat(Hyperrectangle(SA[N(3), N(2)], SA[N(2), N(1)]))
        G3 = @inferred LazySets._genmat_static(Hyperrectangle(SA[N(3), N(2)], SA[N(2), N(1)]))
        @test G1 isa SparseMatrixCSC && G2 isa SMatrix && G3 isa SMatrix
        @test G1 == G2 == G3 == N[2 0; 0 1]
    end

    # high
    res = @inferred high(H)
    @test res isa Vector{N} && res == N[2, 1]
    @test_throws DimensionMismatch high(H, 3)
    res = @inferred high(H, 1)
    @test res isa N && res == N(2)

    # isbounded
    @test @inferred isbounded(H)

    # isboundedtype
    @test @inferred isboundedtype(typeof(H))

    # isconvex
    @test @inferred isconvex(H)

    # isconvextype
    @test @inferred isconvextype(typeof(H))

    # isempty
    @test !(@inferred isempty(H))
    @test_broken @inferred isempty(H, true)  # TODO make this type-stable
    res, w = isempty(H, true)
    @test !res && w isa Vector{N} && w ∈ H

    # isflat
    @test !(@inferred isflat(H))
    @test @inferred isflat(H0)
    @test @inferred isflat(Hyperrectangle(N[1, 2], N[1, 0]))

    # isoperation
    @test !(@inferred isoperation(H))

    # isoperationtype
    @test !(@inferred isoperationtype(typeof(H)))

    # ispolyhedral
    @test @inferred ispolyhedral(H)

    # ispolyhedraltype
    @test @inferred ispolyhedraltype(typeof(H))

    # ispolytopic
    @test @inferred ispolytopic(H)

    # ispolytopictype
    @test @inferred ispolytopictype(typeof(H))

    # isuniversal
    @test !(@inferred isuniversal(H))
    @test_broken @inferred isuniversal(H, true)  # TODO make this type-stable
    res, w = isuniversal(H, true)
    @test !res && w isa Vector{N} && w ∉ H

    # low
    res = @inferred low(H)
    @test res isa Vector{N} && res == N[0, -3]
    @test_throws DimensionMismatch low(H, 3)
    res = @inferred low(H, 1)
    @test res isa N && res == N(0)

    # ngens
    @static if VERSION >= v"1.12"
        @test (@inferred ngens(H)) == 2
    end
    @test ngens(H) == 2
    # degenerate case
    @static if VERSION >= v"1.12"
        @test (@inferred ngens(H)) == 2
    end
    @test ngens(H0) == 0

    # norm
    @test_throws ArgumentError norm(H, N(1 // 2))
    for res in ((@inferred norm(H)), @inferred norm(H, Inf))
        if N <: AbstractFloat
            @test res isa N && res == N(3)
        else
            @test res == 3.0
        end
    end
    res = @inferred norm(H, 2)
    if N <: AbstractFloat
        @test res isa N
    end
    @test res == norm(N[2, -3], 2)

    # radius
    @test_throws ArgumentError radius(H, N(1 // 2))
    for res in ((@inferred radius(H)), @inferred radius(H, Inf))
        if N <: AbstractFloat
            @test res isa N && res == N(2)
        else
            @test res == 2.0
        end
    end
    res = @inferred radius(H, 2)
    if N <: AbstractFloat
        @test res isa N
    end
    @test res == norm(N[1, 2], 2)

    # radius_hyperrectangle
    res = @inferred radius_hyperrectangle(H)
    @test res isa AbstractVector{N} && res == N[1, 2]
    @test_throws DimensionMismatch radius_hyperrectangle(H, 3)
    res = @inferred radius_hyperrectangle(H, 1)
    @test res isa N && res == N(1)

    # rectify
    H2 = @inferred rectify(H)
    @test isidentical(H2, Hyperrectangle(N[1, 1 // 2], N[1, 1 // 2]))

    # reflect
    H2 = @inferred reflect(H)
    @test isidentical(H2, Hyperrectangle(N[-1, 1], r))

    # remove_redundant_generators
    @test isidentical((@inferred remove_redundant_generators(H)), H)

    # singleton_list
    res = @inferred singleton_list(H)
    @test res isa Vector{Singleton{N,Vector{N}}}
    @test ispermutation(res,
                        [Singleton(N[2, 1]), Singleton(N[0, 1]),
                         Singleton(N[0, -3]), Singleton(N[2, -3])])

    # togrep
    Z = @inferred togrep(H)
    @test Z isa Zonotope{N} && isequivalent(Z, H)

    # tosimplehrep
    A, b = @inferred tosimplehrep(H)
    @test A isa Matrix{N} && b isa Vector{N} && size(A) == (4, 2) && length(b) == 4
    # no precise test here; the `tosimplehrep` implementation is tested elsewhere

    # triangulate_faces
    @test_throws DimensionMismatch triangulate_faces(H)
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :GeometryBasics)
        res = @inferred triangulate_faces(H3)
        @test res isa Tuple{Matrix{Float32},Vector{Tuple{Int,Int,Int}}}
        @test length(res[2]) == 12
        # no precise test here; the `triangulate_faces` implementation is tested elsewhere
    end

    # vertices_list
    vlist = @inferred vertices_list(H)
    @test vlist isa Vector{Vector{N}} &&
          ispermutation(vlist, [N[2, 1], N[0, 1], N[0, -3], N[2, -3]])
    # degenerate case
    vlist = @inferred vertices_list(H0)
    @test vlist isa Vector{Vector{N}} && vlist == [N[0, 0]]
    # degenerate case: efficient handling (#92; would not be able to compute all 2^100 vertices)
    c2 = fill(N(1), 100)
    r2 = SingleEntryVector(1, 100, N(1))
    @test ispermutation((@inferred vertices_list(Hyperrectangle(c2, r2))), [c2 - r2, c2 + r2])

    # vertices
    vlist = collect(@inferred vertices(H))
    @test vlist isa Vector{Vector{N}} && ispermutation(vlist, vertices_list(H))

    # volume
    res = @inferred volume(H)
    @test res isa N && res == N(8)
    res = @inferred volume(H0)
    @test res isa N && res == N(0)
    res = @inferred volume(H3)
    @test res isa N && res == N(48)

    # affine_map
    @test_throws DimensionMismatch affine_map(ones(N, 2, 1), H, N[1, 1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 2), H, N[1])
    X = @inferred affine_map(ones(N, 1, 2), H, N[1])
    @test X == Zonotope(N[1], hcat(N(3)))
    X = @inferred affine_map(ones(N, 3, 2), H, N[1, 2, 3])
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test isequivalent(X, Zonotope(N[1, 2, 3], hcat(N[3, 3, 3])))
    end

    # distance (between point and set)
    @test_throws DimensionMismatch distance(H, N[0])
    @test_throws ArgumentError distance(H, N[0, 0]; p=N(1 // 2))
    xs = [N[1, 0], N[5, 1], N[1, 3], N[-3, -1], N[-1, -4]]
    ys = [N[1, 0], N[2, 1], N[1, 1], N[0, -1], N[0, -3]]  # closest points in H
    for (x, y) in zip(xs, ys)
        if N <: AbstractFloat
            @test (@inferred distance(x, H)) == (@inferred distance(H, x)) ==
                  @inferred distance(x, y)
        else
            @test_broken @inferred distance(x, H)  # TODO make this type-stable
            @test distance(x, H) == distance(H, x) == distance(x, y)
        end
        for p in N[1, 2, Inf]
            if N <: AbstractFloat
                @test (@inferred distance(x, H; p=p)) == (@inferred distance(H, x; p=p)) ==
                      @inferred distance(x, y; p=p)
            else
                @test distance(x, H; p=p) == distance(H, x; p=p) == distance(x, y; p=p)
            end
        end
    end

    # exponential_map
    @test_throws DimensionMismatch exponential_map(ones(N, 1, 1), H)
    @test_throws DimensionMismatch exponential_map(ones(N, 2, 1), H)
    X = @inferred exponential_map(N[1 0; 0 2], H)
    @test X isa Zonotope && center(X) == N[exp(1), -exp(2)] && genmat(X) == N[exp(1) 0; 0 2*exp(2)]

    # in
    @test_throws DimensionMismatch N[0] ∈ H
    @test (@inferred N[1, -2] ∈ H) && (@inferred N[2, -3] ∈ H) && (@inferred N[3, 2] ∉ H)
    @test (@inferred N[0, 0] ∈ H0) && (@inferred N[2, 3] ∉ H0)
    # robust membership (#1576)
    H2 = Hyperrectangle(N[1.68, 0.73, 0.64], N[0.46, 0.24, 1.38])
    @test @inferred σ(ones(N, 3), H2) ∈ H2

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0], H)
    @test_throws ArgumentError is_interior_point(N[0, 0], H; ε=N(0))
    @test_throws ArgumentError is_interior_point(N[0, 0], H; p=N(1 // 2))
    if N <: AbstractFloat
        @test @inferred is_interior_point(N[1, 0], H)
        @test !(@inferred is_interior_point(N[5, 8], H))
        @test !(@inferred is_interior_point(N[2, 1], H))
    else
        @test_throws ArgumentError is_interior_point(N[1, 0], H)
        @test is_interior_point(N[1, 0], H; ε=1 // 100)
        @test !(@inferred is_interior_point(N[5, 8], H; ε=1 // 100))
        @test !(@inferred is_interior_point(N[2, 1], H; ε=1 // 100))
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([0.0, 0.0], H)
    end

    # linear_map
    @test_throws DimensionMismatch linear_map(ones(N, 1, 1), H)
    X = @inferred linear_map(N[1 0; 0 2], H)
    @test X isa Zonotope{N} && X == Zonotope(N[1, -2], N[1 0; 0 4])
    X = @inferred linear_map(zeros(N, 1, 2), H)  # zero map
    @test X isa Zonotope{N} && X == Zonotope(N[0], zeros(N, 1, 0))

    # linear_map_inverse
    @test_throws AssertionError LazySets.linear_map_inverse(ones(N, 1, 1), H)
    # invertible map
    X = LazySets.linear_map_inverse(N[1 0; 0 1], H)
    @test X isa Zonotope{N} && X == Zonotope(N[1, -1], N[1 0; 0 2])
    # noninvertible map
    M = ones(N, 2, 1)
    X = LazySets.linear_map_inverse(M, H)
    A, b = tosimplehrep(H)
    Y = HPolytope(A * M, b)
    @test X isa LazySet{N} && isequivalent(X, Y)

    # permute
    @test_throws DimensionMismatch permute(H, [1])
    @test_throws DimensionMismatch permute(H, [1, -1])
    @test_throws DimensionMismatch permute(H, [1, 3])
    @test_throws ArgumentError permute(H, [1, 1])
    @test isidentical((@inferred permute(H, [1, 2])), H)
    H2 = @inferred permute(H, [2, 1])
    @test isidentical(H2, Hyperrectangle(N[-1, 1], N[2, 1]))

    # project
    @test_throws DimensionMismatch project(H, [1, 2, 3])
    @test_throws DimensionMismatch project(H, [1, -1])
    @test_throws DimensionMismatch project(H, [1, 3])
    @test_throws ArgumentError project(H, [1, 1])
    H2 = @inferred project(H, [1])
    @test isidentical(H2, Hyperrectangle(N[1], N[1]))
    H2 = @inferred project(Hyperrectangle(N[4, 3, 2, 1], N[8, 7, 6, 5]), [2, 4])
    @test isidentical(H2, Hyperrectangle(N[3, 1], N[7, 5]))

    # sample
    res = @inferred sample(H)
    @test res isa Vector{N} && res ∈ H
    res = @inferred sample(H, 2)
    @test res isa Vector{Vector{N}} && length(res) == 2 && all(x ∈ H for x in res)

    # scale
    H2 = @inferred scale(N(2), H)
    @test isidentical(H2, Hyperrectangle(N[2, -2], N[2, 4]))
    # degenerate case
    H2 = @inferred scale(N(0), H)
    @test isidentical(H2, Hyperrectangle(N[0, 0], N[0, 0]))
    # scale!
    H2 = copy(H)
    @inferred scale!(N(2), H2)
    @test isidentical(H2, Hyperrectangle(N[2, -2], N[2, 4]))
    # degenerate case
    H2 = copy(H)
    @inferred scale!(N(0), H2)
    @test isidentical(H2, Hyperrectangle(N[0, 0], N[0, 0]))
    # negative factor
    H2 = copy(H)
    @inferred scale!(N(-2), H2)
    @test isidentical(scale(N(-2), H), H2) && isidentical(H2, Hyperrectangle(N[-2, 2], N[2, 4]))

    # split
    @test_throws ArgumentError split(H, [0, 2])
    S = [Hyperrectangle(N[1 // 2, -2], N[1 // 2, 1]),
         Hyperrectangle(N[3 // 2, -2], N[1 // 2, 1]),
         Hyperrectangle(N[1 // 2, 0], N[1 // 2, 1]),
         Hyperrectangle(N[3 // 2, 0], N[1 // 2, 1])]
    @test ispermutation((@inferred split(H, [2, 2])), S)
    @static if isdefined(@__MODULE__, :StaticArrays)
        # static vectors
        using StaticArrays: SA

        H2 = Hyperrectangle(SA[N(1), N(-1)], SA[N(1), N(2)])
        S2 = @inferred split(H2, [2, 2])
        @test S2 isa Vector{typeof(H2)} && ispermutation(S2, S)
    end
    # split at given bounds
    @test ispermutation(split(H, [N[1], N[-1]]), S)
    S2 = @inferred split(H, [N[1, 3 // 2], N[-1, 0]])
    @test ispermutation(S2,
                        [Hyperrectangle(; low=N[0, -3], high=N[1, -1]),
                         Hyperrectangle(; low=N[0, -1], high=N[1, 0]),
                         Hyperrectangle(; low=N[0, 0], high=N[1, 1]),
                         Hyperrectangle(; low=N[1, -3], high=N[3 // 2, -1]),
                         Hyperrectangle(; low=N[1, -1], high=N[3 // 2, 0]),
                         Hyperrectangle(; low=N[1, 0], high=N[3 // 2, 1]),
                         Hyperrectangle(; low=N[3 // 2, -3], high=N[2, -1]),
                         Hyperrectangle(; low=N[3 // 2, -1], high=N[2, 0]),
                         Hyperrectangle(; low=N[3 // 2, 0], high=N[2, 1])])

    # support_function
    @test_throws DimensionMismatch ρ(N[1], H)
    res = @inferred ρ(N[1, 1], H)
    @test res isa N && res == N(3)
    res = @inferred ρ(N[-1, -1], H)
    @test res isa N && res == N(3)
    # SingleEntryVector
    res = @inferred ρ(SingleEntryVector(2, 2, N(-1)), H)
    @test res isa N && res == N(3)

    # support_vector
    @test_throws DimensionMismatch σ(N[1], H)
    res = @inferred σ(N[1, 1], H)
    @test res isa Vector{N} && res == N[2, 1]
    res = @inferred σ(N[-1, -1], H)
    @test res isa Vector{N} && res == N[0, -3]

    # translate
    @test_throws DimensionMismatch translate(H, N[1])
    H2 = @inferred translate(H, N[1, 2])
    @test isidentical(H2, Hyperrectangle(N[2, 1], N[1, 2]))
    # translate!
    @test_throws DimensionMismatch translate!(H, N[1])
    H2 = copy(H)
    @inferred translate!(H2, N[1, 2])
    @test isidentical(H2, Hyperrectangle(N[2, 1], N[1, 2]))

    # cartesian_product
    H2a = Hyperrectangle(N[3], N[4])
    H2 = @inferred cartesian_product(H1, H2a)
    @test isidentical(H2, Hyperrectangle(N[0, 3], N[1, 4]))
    H2 = @inferred cartesian_product(H2a, H1)
    @test isidentical(H2, Hyperrectangle(N[3, 0], N[4, 1]))

    # convex_hull (binary)
    @test_throws DimensionMismatch convex_hull(H, H3)
    @test_broken @inferred convex_hull(H, H)  # TODO make this type-stable
    X = convex_hull(H, H)
    @test X isa LazySet{N} && isequivalent(X, H)
    H2 = Hyperrectangle(N[2, 0], N[1, 2])
    Y = VPolygon([N[3, 2], N[1, 2], N[0, 1], N[0, -3], N[2, -3], N[3, -2]])
    for X in (convex_hull(H, H2), convex_hull(H2, H))
        @test X isa LazySet{N} && isequivalent(X, Y)
    end

    # difference
    @test_throws DimensionMismatch difference(H, H3)
    # full set
    @static if isdefined(@__MODULE__, :IntervalBoxes)
        X = difference(H, Hyperrectangle(N[4, 3], N[1, 1]))
        @test X isa LazySet{N} && isequivalent(X, H)
        # empty set
        X = difference(H, H)
        @test isequivalent(X, EmptySet{N}(2))
        # one hyperrectangle
        X = difference(H, Hyperrectangle(N[2, -1], N[1, 2]))
        @test X isa LazySet{N} && isequivalent(X, Hyperrectangle(N[1 // 2, -1], N[1 // 2, 2]))
        # multiple hyperrectangles (result is not unique; reylying on implementation here)
        X = difference(H, Hyperrectangle(N[2, 0], N[1, 2]))
        @test X isa UnionSetArray{N} &&
              array(X) == [Hyperrectangle(N[1 // 2, -1], N[1 // 2, 2]),
                           Hyperrectangle(N[3 // 2, -5 // 2], N[1 // 2, 1 // 2])]
    end

    # distance (between two sets)
    @test_throws DimensionMismatch distance(H, H3)
    @test_throws ArgumentError distance(H, H; p=N(1 // 2))
    res = @inferred distance(H, H)
    if N <: AbstractFloat
        @test res isa N && res == N(0)
    else
        @test res isa Float64 && res == 0.0
    end
    H2 = Hyperrectangle(N[4, 3], N[1, 1])
    x = N[2, 1]
    y = N[3, 2]
    @test (@inferred distance(H, H2)) == (@inferred distance(H2, H)) == @inferred distance(x, y)
    for p in N[1, 2, Inf]
        @test (@inferred distance(H, H2; p=p)) == (@inferred distance(H2, H; p=p)) ==
              @inferred distance(x, y; p=p)
    end

    # exact_sum
    @test_throws DimensionMismatch exact_sum(H, H3)
    H2 = Hyperrectangle(N[4, -2], N[2, 4])
    for H2b in ((@inferred exact_sum(H, H2)), @inferred exact_sum(H2, H))
        @test H2b isa Hyperrectangle{N} &&
              isequivalent(H2b, Hyperrectangle(N[5, -3], N[3, 6]))
    end

    # intersection
    @test_throws DimensionMismatch intersection(H, H3)
    # disjoint
    X = intersection(H, Hyperrectangle(N[4, 0], N[1, 1]))
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    # overlapping
    H2 = intersection(H, Hyperrectangle(N[2, 0], N[1, 2]))
    @test isidentical(H2, Hyperrectangle(N[3 // 2, -1 // 2], N[1 // 2, 3 // 2]))

    # isapprox
    @test @inferred H ≈ H
    res = (H ≈ translate(H, N[1 // 100000000, 0]))
    if N <: AbstractFloat
        @test res  # below default tolerance for AbstractFloat
    else
        @test !res  # zero default tolerance for Rational
    end
    @test !(@inferred H ≈ translate(H, N[1 // 1000, 0]))  # above default tolerance for all types
    @test !(@inferred H ≈ H3) && !(@inferred H3 ≈ H) && !(@inferred H ≈ P) && !(@inferred P ≈ H)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(H, H3)
    # disjoint
    H2 = Hyperrectangle(N[4, 3], N[1, 1])
    for (H2a, H2b) in ((H, H2), (H2, H))
        @test_broken @inferred isdisjoint(H2a, H2b)  # TODO make this type-stable
        @test isdisjoint(H2a, H2b)
        @test_broken @inferred isdisjoint(H2a, H2b, true)  # TODO make this type-stable
        res, w = isdisjoint(H2a, H2b, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    # overlapping
    H2 = Hyperrectangle(N[2, 2], N[1, 2])
    for (H2a, H2b) in ((H, H), (H, H2), (H2, H))
        @test !isdisjoint(H2a, H2b)
        res, w = isdisjoint(H2a, H2b, true)
        @test !res && w isa Vector{N} && w ∈ H2a && w ∈ H2b
    end
    # flat
    H2 = Hyperrectangle(N[2, 3], zeros(N, 2))
    for (H2a, H2b) in ((H0, H2), (H2, H0))
        @test isdisjoint(H2a, H2b)
        res, w = isdisjoint(H2a, H2b, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    @test !isdisjoint(H0, H0)
    res, w = isdisjoint(H0, H0, true)
    @test !res && w isa Vector{N} && w ∈ H0
    # tolerance
    if N == Float64
        H2 = Hyperrectangle(N[1 + 1e-9, -1], N[1, 2])
        @test !isdisjoint(H, H2)
        r = _rtol(N)
        @assert r > N(1e-10) "default tolerance changed; adapt test"
        set_rtol(N, N(1e-10))
        @test_broken isdisjoint(H, H2)  # cannot adapt tolerance of LP solver
        # restore tolerance
        set_rtol(N, r)
    end

    # isequal
    @test @inferred H == H
    @test (@inferred H != H3) && (@inferred H3 != H) && (@inferred H != P) && @inferred P != H

    # isequivalent
    @test_throws DimensionMismatch isequivalent(H, H3)
    @test_broken @inferred isequivalent(H, H)  # TODO make this type-stable
    @test isequivalent(H, H)
    @test !isequivalent(H, Hyperrectangle(N[1, 2], N[1, 1]))
    @test isequivalent(H, P) && isequivalent(P, H)

    # isstrictsubset
    @test_throws DimensionMismatch H ⊂ H3
    for X in (H, P)
        @test_broken !(@inferred X ⊂ H)  # TODO make this type-stable
        @test !(X ⊂ H)
        @test_broken @inferred ⊂(X, H, true)  # TODO make this type-stable
        res, w = ⊂(X, H, true)
        @test !res && w isa Vector{N} && isempty(w)
    end
    H2 = Hyperrectangle(N[1, 2], N[1, 1])
    @test !(H2 ⊂ H)
    res, w = ⊂(H2, H, true)
    @test !res && w isa Vector{N} && w ∈ H2 && w ∉ H
    H2 = Hyperrectangle(N[1, -1], N[2, 3])
    @test H ⊂ H2
    res, w = ⊂(H, H2, true)
    @test res && w isa Vector{N} && w ∉ H && w ∈ H2

    # issubset
    @test_throws DimensionMismatch H ⊆ H3
    for X in (H, P, Hyperrectangle(N[1, -1], N[2, 3]))
        @test_broken @inferred H ⊆ X  # TODO make this type-stable
        @test H ⊆ X
        @test_broken @inferred ⊆(H, X, true)  # TODO make this type-stable
        res, w = ⊆(H, X, true)
        @test res && w isa Vector{N} && w == N[]
    end
    H2 = Hyperrectangle(N[3, 0], N[1, 1])
    @test H ⊈ H2
    res, w = ⊆(H, H2, true)
    @test !res && w isa Vector{N} && w ∈ H && w ∉ H2

    # linear_combination
    @test_throws DimensionMismatch linear_combination(H, H3)
    @test_broken linear_combination(H, Xnc) isa LazySet{N}  # TODO implement `linear_combination` for non-convex sets
    @test_broken linear_combination(Xnc, H) isa LazySet{N}
    @test_broken @inferred linear_combination(H, H)  # TODO make this type-stable
    for X in (linear_combination(H, H), linear_combination(H, P), linear_combination(P, H))
        @test X isa LazySet{N} && isequivalent(X, H)
    end
    H2 = Hyperrectangle(N[2, 0], N[1, 2])
    Y = VPolygon([N[3, 2], N[1, 2], N[0, 1], N[0, -3], N[2, -3], N[3, -2]])
    for X in (linear_combination(H, H2), linear_combination(H2, H))
        @test X isa LazySet{N} && isequivalent(X, Y)
    end

    # minkowski_difference
    @test_throws DimensionMismatch minkowski_difference(H, H3)
    # empty difference
    H2 = Hyperrectangle(N[1, 2], N[4, 6])
    X = minkowski_difference(H, H2)
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    # equivalent sets: only the origin remains
    H2 = minkowski_difference(H, H)
    @test isidentical(H2, Hyperrectangle(N[0, 0], zeros(N, 2)))
    for X in (minkowski_difference(H, P), minkowski_difference(P, H))
        vlist = vertices_list(X)
        @test length(vlist) == 1 && iszero(vlist)
    end
    # nonempty difference
    H2 = minkowski_difference(H, Hyperrectangle(N[1, 0], N[1, 1]))
    @test isidentical(H2, Hyperrectangle(N[0, -1], N[0, 1]))

    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(H, H3)
    H2 = Hyperrectangle(N[4, -2], N[2, 4])
    for H2b in ((@inferred minkowski_sum(H, H2)), @inferred minkowski_sum(H2, H))
        @test H2b isa Hyperrectangle{N} &&
              isequivalent(H2b, Hyperrectangle(N[5, -3], N[3, 6]))
    end
end

for N in @tN([Float64, Float32])
    H = Hyperrectangle(N[1, -1], N[1, 2])

    # rand
    H2 = rand(Hyperrectangle; N=N)
    @test H2 isa Hyperrectangle{N} && dim(H2) == 2
    H2 = rand(Hyperrectangle; N=N, dim=3)
    @test H2 isa Hyperrectangle{N} && dim(H2) == 3

    # rationalize
    @test_broken @inferred rationalize(H)  # TODO make this type-stable
    H2 = rationalize(H)
    T = Rational{Int64}
    @test H2 isa Hyperrectangle{T} && isidentical(H2, Hyperrectangle(T[1, -1], T[1, 2]))
end

for N in @tN([Float64, Rational{Int}])
    H = Hyperrectangle(N[1, -1], N[1, 2])

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        X = polyhedron(H)
        @test X isa Polyhedra.DefaultPolyhedron{N}
        @test ispermutation(constraints_list(convert(HPolytope, X)), constraints_list(H))
    end

    # triangulate
    @static if isdefined(@__MODULE__, :MiniQhull)
        X = @inferred triangulate(H)
        @test X isa UnionSetArray{N,<:VPolytope{N}} && length(X) == 2
        # result is not unique; reylying on implementation here
        for Xi in X
            @test ispermutation(vertices_list(Xi), [N[2, 1], N[0, -3], N[2, -3]]) ||
                  ispermutation(vertices_list(Xi), [N[2, 1], N[0, -3], N[0, 1]])
        end
    end
end
