using LazySets, Test
using LazySets.ReachabilityBase.Arrays: ispermutation, SingleEntryVector
using LazySets.ReachabilityBase.Comparison: set_rtol, _rtol, _ztol
IA = LazySets.IA
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

function isidentical(::Interval, ::Interval)
    return false
end

function isidentical(X1::Interval{N}, X2::Interval{N}) where {N}
    return X1 == X2
end

for N in @tN([Float64, Float32, Rational{Int}])
    # auxiliary sets
    S2 = Singleton(N[0, 0])  # 2D set
    B = BallInf(N[1], N(1))  # equivalent set
    Xnc = UnionSet(B, BallInf(N[4], N(1)))  # nonconvex set

    # default constructor from IntervalArithmetic.Interval
    itv = IA.interval(N(0), N(2))
    X = @inferred Interval(itv)
    @test X isa Interval{N} && IA.isequal_interval(X.dat, itv)

    # constructors from two numbers, from a vector, and with promotion
    for Y in ((@inferred Interval(N(0), N(2))), (@inferred Interval(N[0, 2])),
              @inferred Interval(0, N(2)))
        @test isidentical(Y, X)
    end

    # constructor from a number
    X0 = @inferred Interval(N(0))
    @test X0 isa Interval{N} && X0 == Interval(N(0), N(0))

    # unbounded intervals are caught
    @test_throws AssertionError Interval(-Inf, zero(N))
    @test_throws AssertionError Interval(zero(N), Inf)

    # convert
    # to and from IntervalArithmetic.Interval
    Y = @inferred convert(IA.Interval, X)
    @test IA.isequal_interval(Y, X.dat)
    Z = @inferred convert(Interval, Y)
    @test isidentical(Z, X)
    # from hyperrectangular set
    @test_throws AssertionError convert(Interval, Hyperrectangle(N[0, 0], N[1, 1]))
    Y = Hyperrectangle(N[1], N[1])
    Z = @inferred convert(Interval, Y)
    @test isidentical(Z, X)
    # from Rectification of Interval
    Y = Interval(N(-1), N(2))
    Z = @inferred convert(Interval, Rectification(Y))
    @test isidentical(Z, X)
    # from MinkowskiSum of Intervals
    Z = @inferred convert(Interval, Y + Y)
    @test isidentical(Z, Interval(N(-2), N(4)))
    # from LazySet
    M = hcat(N[1])
    Y = @inferred convert(Interval, M * X)
    @test isidentical(Y, X)
    # from non-convex set
    Y = UnionSet(Interval(1, 2), Interval(3, 4))
    @test_throws AssertionError convert(Interval, Y)

    # an_element
    x = @inferred an_element(X)
    @test x isa Vector{N} && length(x) == 1 && N(0) <= x[1] <= N(2)

    # area
    @test_throws DimensionMismatch area(X)

    # center
    c = @inferred center(X)
    @test c isa Vector{N} && c == [N(1)]
    v = @inferred center(X, 1)
    @test v isa N && v == N(1)
    @test_throws DimensionMismatch center(X, 2)

    # chebyshev_center_radius
    c, r = @inferred chebyshev_center_radius(X)
    @test c isa Vector{N} && c == [N(1)]
    @test r isa N && r == N(1)

    # complement
    Y = @inferred complement(X)
    @test Y isa UnionSet{N}
    @test ispermutation(array(Y), [HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-2))])

    # concretize
    Y = @inferred concretize(X)
    @test isidentical(Y, X)

    # constrained_dimensions
    @test (@inferred constrained_dimensions(X)) == 1:1

    # constraints_list
    cs = @inferred constraints_list(X)
    @test ispermutation(cs, [HalfSpace(N[1], N(2)), HalfSpace(N[-1], N(0))])

    # constraints
    cs2 = collect(@inferred constraints(X))
    @test ispermutation(cs2, cs)

    # convex_hull (unary)
    Y = @inferred convex_hull(X)
    @test isidentical(Y, X)

    # copy
    Y = @inferred copy(X)
    @test isidentical(Y, X)

    # diameter
    @test_throws ArgumentError diameter(X, N(1 // 2))
    for res in ((@inferred diameter(X)), (@inferred diameter(X, Inf)), @inferred diameter(X, 2))
        @test res isa N && res == N(2)
    end

    # dim
    @test (@inferred dim(X)) == 1

    # eltype
    @test (@inferred eltype(X)) == N
    @test (@inferred eltype(typeof(X))) == N

    # extrema
    res = @inferred extrema(X)
    @test res isa Tuple{Vector{N},Vector{N}} && res[1] == N[0] && res[2] == N[2]
    @test_throws DimensionMismatch extrema(X, 2)
    res = @inferred extrema(X, 1)
    @test res isa Tuple{N,N} && res[1] == N(0) && res[2] == N(2)

    # generators
    @test collect(@inferred generators(X)) == [N[1]]
    # degenerate case
    gens = collect(@inferred generators(X0))
    @test gens isa Vector{SingleEntryVector{N}} && isempty(gens)

    # genmat
    @test (@inferred genmat(X)) == hcat(N[1])
    # degenerate case
    gens = @inferred genmat(X0)
    @test gens isa Matrix{N} && gens == Matrix{N}(undef, 1, 0)

    # high
    res = @inferred high(X)
    @test res isa Vector{N} && res == N[2]
    @test_throws DimensionMismatch high(X, 2)
    res = @inferred high(X, 1)
    @test res isa N && res == N(2)

    # isbounded
    @test @inferred isbounded(X)

    # isboundedtype
    @test @inferred isboundedtype(typeof(X))

    # isconvex
    @test @inferred isconvex(X)

    # isconvextype
    @test @inferred isconvextype(typeof(X))

    # isempty
    @test !(@inferred isempty(X))
    @test_broken @inferred isempty(X, true)  # TODO make this type-stable
    res, w = isempty(X, true)
    @test !res && w ∈ X && w isa Vector{N}

    # isflat
    ztol = _ztol(N)  # default absolute zero tolerance
    @test @inferred isflat(Interval(N(0), ztol))
    @test !(@inferred isflat(Interval(N(0), 2 * ztol + N(1 // 100))))

    # isoperation
    @test !(@inferred isoperation(X))

    # isoperationtype
    @test !(@inferred isoperationtype(typeof(X)))

    # ispolyhedral
    @test @inferred ispolyhedral(X)

    # ispolyhedraltype
    @test @inferred ispolyhedraltype(typeof(X))

    # ispolytopic
    @test @inferred ispolytopic(X)

    # ispolytopictype
    @test @inferred ispolytopictype(typeof(X))

    # isuniversal
    @test !(@inferred isuniversal(X))
    @test_broken @inferred isuniversal(X, true)  # TODO make this type-stable
    res, w = isuniversal(X, true)
    @test !res && w isa Vector{N} && w ∉ X

    # low
    res = @inferred low(X)
    @test res isa Vector{N} && res == N[0]
    @test_throws DimensionMismatch low(X, 2)
    res = @inferred low(X, 1)
    @test res isa N && res == N(0)

    # min
    v = @inferred min(X)
    @test v isa N && v == N(0)

    # max
    v = @inferred max(X)
    @test v isa N && v == N(2)

    # ngens
    @test (@inferred ngens(X)) == 1
    # degenerate case
    @test (@inferred ngens(X0)) == 0

    # norm
    @test_throws ArgumentError norm(X, N(1 // 2))
    for Y in (X, Interval(N(-2), N(1)))
        for res in ((@inferred norm(Y)), (@inferred norm(Y, Inf)), @inferred norm(Y, 2))
            @test res isa N && res == N(2)
        end
    end

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        Y = @inferred polyhedron(X)
        @test Y isa Polyhedra.Interval{N} && ispermutation(Y.vrep.points.points, [[N(0)], [N(2)]])
    end

    # radius
    @test_throws ArgumentError radius(X, N(1 // 2))
    for res in ((@inferred radius(X)), (@inferred radius(X, Inf)), @inferred radius(X, 2))
        @test res isa N && res == N(1)
    end

    # radius_hyperrectangle
    r = @inferred radius_hyperrectangle(X)
    @test r isa Vector{N} && r == [N(1)]
    v = @inferred radius_hyperrectangle(X, 1)
    @test v isa N && v == N(1)

    # rectify
    Y = Interval(N(-1), N(2))
    Z = @inferred rectify(Y)
    @test isidentical(Z, X)
    Y = Interval(N(-2), N(-1))
    Z = @inferred rectify(Y)
    @test isidentical(Z, X0)
    Y = Interval(N(1), N(2))
    Z = @inferred rectify(Y)
    @test isidentical(Y, Z)

    # reflect
    Y = @inferred reflect(X)
    @test isidentical(Y, Interval(N(-2), N(0)))

    # singleton_list
    res = @inferred singleton_list(X)
    @test res isa Vector{Singleton{N,Vector{N}}}
    @test ispermutation(res, [Singleton(N[0]), Singleton(N[2])])

    # togrep
    Z = @inferred togrep(X)
    @test Z isa Zonotope{N} && isequivalent(Z, X)

    # tosimplehrep
    A, b = @inferred tosimplehrep(X)
    @test A isa Matrix{N} && b isa Vector{N}
    @test (A == hcat(N[1, -1]) && b == N[2, 0]) || (A == hcat(N[-1, 1]) && b == N[0, 2])

    # triangulate
    @static if isdefined(@__MODULE__, :MiniQhull)
        Y = @inferred triangulate(X)  # TODO does it make sense to triangulate a set with two vertices?
        @test Y isa UnionSetArray{N} && length(Y) == 1
        Y = Y[1]
        @test Y isa VPolytope{N}
        res = vertices_list(Y)
        @test ispermutation(res, [N[0], N[2]])
    end

    # triangulate_faces
    @test_throws DimensionMismatch triangulate_faces(X)

    # vertices_list
    res = @inferred vertices_list(X)
    @test res isa Vector{Vector{N}} && ispermutation(res, [N[0], N[2]])
    # degenerate case
    res = @inferred vertices_list(X0)
    @test res isa Vector{Vector{N}} && res == [[N(0)]]

    # vertices
    res = collect(@inferred vertices(X))
    @test res isa Vector{Vector{N}} && ispermutation(res, vertices_list(X))

    # volume
    @test (@inferred volume(X)) == N(2)

    # affine_map
    @test_throws DimensionMismatch affine_map(ones(N, 1, 2), X, N[1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 1), X, N[1])
    Y = affine_map(ones(N, 1, 1), X, N[1])
    @test isidentical(Y, Interval(N(1), N(3)))
    Y = affine_map(ones(N, 2, 1), X, N[1, 2])
    @test Y isa LazySet{N} && isequivalent(Y, LineSegment(N[1, 2], N[3, 4]))

    # distance (between point and set)
    @test_throws DimensionMismatch distance(X, N[0, 0])
    @test_throws ArgumentError distance(X, N[0]; p=N(1 // 2))
    for (x, v) in ((N[1], N(0)), (N[4], N(2)))
        for res in (distance(X, x), distance(x, X))
            @test res == v
            if N <: AbstractFloat
                @inferred distance(X, x)
                @inferred distance(x, X)
                @test res isa N
            else
                @test_broken @inferred distance(X, x)  # TODO make this type-stable
            end
        end
    end

    # exponential_map
    @test_throws DimensionMismatch exponential_map(ones(N, 2, 2), X)
    @test_throws DimensionMismatch exponential_map(ones(N, 2, 1), X)

    # in
    @test_throws DimensionMismatch N[0, 0] ∈ X
    @test (@inferred N[1] ∈ X) && N[2] ∈ X && N[3] ∉ X
    # number in interval is invalid
    @test_throws MethodError N(1) ∈ X

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0, 0], X)
    @test_throws ArgumentError is_interior_point(N[0], X; ε=N(0))
    @test_throws ArgumentError is_interior_point(N[0], X; p=N(1 // 2))
    if N <: AbstractFloat
        @test @inferred is_interior_point(N[1], X)
        @test !(@inferred is_interior_point(N[2], X))
        @test !(@inferred is_interior_point(N[3], X))
        @test @inferred is_interior_point(N[1], X; p=N(2))
        @test @inferred is_interior_point(N[1], X; p=N(2), ε=N(1 // 2))
        @test !(@inferred is_interior_point(N[1], X; ε=N(1)))
    else
        @test_throws ArgumentError is_interior_point(N[1], X)
        @test @inferred is_interior_point(N[1], X; ε=1 // 100)
        @test !(@inferred is_interior_point(N[2], X; ε=1 // 100))
        @test !(@inferred is_interior_point(N[3], X; ε=1 // 100))
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([0.0], X)
    end

    # linear_map
    @test_throws DimensionMismatch linear_map(ones(N, 2, 2), X)
    Y = linear_map(2 * ones(N, 1, 1), X)
    @test isidentical(Y, Interval(N(0), N(4)))
    Y = linear_map(zeros(N, 1, 1), X)  # zero map
    @test Y isa Interval{N}
    @test isequivalent(Y, Interval(N(0), N(0)))
    Y = linear_map(ones(N, 2, 1), X)  # higher dimension
    @test Y isa LazySet{N} && isequivalent(Y, LineSegment(N[0, 0], N[2, 2]))
    Y = linear_map(zeros(N, 2, 1), X)  # zero map in higher dimension
    @test Y isa LazySet{N} && isequivalent(Y, ZeroSet{N}(2))

    # linear_map_inverse
    @test_throws AssertionError LazySets.linear_map_inverse(ones(N, 2, 2), X)
    Y = LazySets.linear_map_inverse(ones(N, 1, 1), X)
    @test Y isa LazySet{N} && isequivalent(Y, X)
    Y = LazySets.linear_map_inverse(ones(N, 1, 2), X)
    Z = HPolyhedron([HalfSpace(N[1, 1], N(2)), HalfSpace(N[-1, -1], N(0))])
    @test Y isa LazySet{N} && isequivalent(Y, Z)

    # permute
    @test_throws DimensionMismatch permute(X, [1, 1])
    @test_throws DimensionMismatch permute(X, [-1])
    @test_throws DimensionMismatch permute(X, [2])
    Y = @inferred permute(X, [1])
    @test isidentical(Y, X)

    # project
    @test_throws DimensionMismatch project(X, [1, 1])
    @test_throws DimensionMismatch project(X, [-1])
    @test_throws DimensionMismatch project(X, [2])
    Y = @inferred project(X, [1])
    @test isidentical(Y, X)

    # sample
    res = @inferred sample(X)
    @test res isa Vector{N} && res ∈ X
    res = @inferred sample(X, 2)
    @test res isa Vector{Vector{N}} && length(res) == 2 && all(x ∈ X for x in res)

    # scale
    Y = @inferred scale(N(2), X)
    @test isidentical(Y, Interval(N(0), N(4)))
    # degenerate case
    Y = @inferred scale(N(0), X)
    @test isidentical(Y, X0)
    # scale!
    @test_throws MethodError scale!(N(2), X)

    # split
    @test_throws AssertionError split(X, 0)
    @test_throws AssertionError split(X, [4, 4])
    Xs = [Interval(N(0), N(1 // 2)), Interval(N(1 // 2), N(1)),
          Interval(N(1), N(3 // 2)), Interval(N(3 // 2), N(2))]
    Ys = @inferred split(X, 4)
    @test Ys isa Vector{Interval{N}} && Ys == split(X, [4]) == Xs

    # support_function
    @test_throws DimensionMismatch ρ(N[1, 1], X)
    res = @inferred ρ(N[2], X)
    @test res isa N && res == N(4)
    res = @inferred ρ(N[-2], X)
    @test res isa N && res == N(0)

    # support_vector
    @test_throws DimensionMismatch σ(N[1, 1], X)
    res = @inferred σ(N[2], X)
    @test res isa Vector{N} && res == N[2]
    res = @inferred σ(N[-2], X)
    @test res isa Vector{N} && res == N[0]

    # translate
    @test_throws DimensionMismatch translate(X, N[1, 1])
    Y = @inferred translate(X, N[1])
    @test isidentical(Y, Interval(N(1), N(3)))
    # translate!
    @test_throws MethodError translate!(X, N[1, 1])  # TODO this should maybe change
    @test_throws MethodError translate!(X, N[1])  # TODO this should maybe change

    # cartesian_product
    Y = Interval(N(-3), N(1))
    Z = @inferred cartesian_product(X, Y)
    @test Z isa Hyperrectangle{N} && Z == Hyperrectangle(N[1, -1], N[1, 2])
    Z = @inferred cartesian_product(Y, X)
    @test Z isa Hyperrectangle{N} && Z == Hyperrectangle(N[-1, 1], N[2, 1])

    # convex_hull (binary)
    @test_throws DimensionMismatch convex_hull(X, S2)
    @test_throws DimensionMismatch convex_hull(S2, X)
    Y = @inferred convex_hull(X, X)
    @test isidentical(Y, X)
    Y = Interval(N(-3), N(-1))
    Z = Interval(N(-3), N(2))
    for W in ((@inferred convex_hull(X, Y)), @inferred convex_hull(Y, X))
        @test W isa LazySet{N} && isidentical(W, Z)
    end

    # difference
    @test_throws DimensionMismatch difference(X, S2)
    @test_throws DimensionMismatch difference(S2, X)
    # disjoint
    @test isidentical(difference(X, Interval(N(3), N(4))), X)
    # overlapping
    @test isidentical(difference(X, Interval(N(0), N(1))), Interval(N(1), N(2)))
    @test isidentical(difference(X, Interval(N(1), N(2))), Interval(N(0), N(1)))
    # fully covered
    Y = difference(X, X)
    @test Y isa EmptySet{N} && Y == EmptySet{N}(1)
    # cut out
    Y = difference(X, Interval(N(1 // 2), N(1)))
    @test Y isa UnionSet{N}
    @test ispermutation(array(Y), [Interval(N(0), N(1 // 2)), Interval(N(1), N(2))])
    # cut out at a single point
    Y = difference(X, Interval(N(1), N(1)))
    @test Y isa UnionSet{N}
    @test ispermutation(array(Y), [Interval(N(0), N(1)), Interval(N(1), N(2))])

    # distance (between two sets)
    @test_throws DimensionMismatch distance(X, S2)
    @test_throws DimensionMismatch distance(S2, X)
    @test_throws ArgumentError distance(X, X; p=N(1 // 2))
    for (Y, v) in ((Interval(N(-1), N(1)), N(0)), (Interval(N(4), N(5)), N(2)))
        for res in ((@inferred distance(X, Y)), @inferred distance(Y, X))
            @test res isa N && res == v
        end
    end

    # exact_sum
    @test_throws DimensionMismatch exact_sum(X, S2)
    @test_throws DimensionMismatch exact_sum(S2, X)
    Y = Interval(N(3), N(4))
    for Z in ((@inferred exact_sum(X, Y)), @inferred exact_sum(Y, X))
        @test isidentical(Z, Interval(N(3), N(6)))
    end

    # intersection
    @test_throws DimensionMismatch intersection(X, S2)
    @test_throws DimensionMismatch intersection(S2, X)
    # disjoint
    Y = intersection(X, Interval(N(3), N(4)))
    @test Y isa EmptySet{N} && Y == EmptySet{N}(1)
    # overlapping
    Y = intersection(X, Interval(N(1), N(3)))
    @test isidentical(Y, Interval(N(1), N(2)))

    # isapprox
    @test @inferred X ≈ X
    res = (@inferred X ≈ translate(X, N[1 // 100000000]))
    if N <: AbstractFloat
        @test res  # below default tolerance for AbstractFloat
    else
        @test !res  # zero default tolerance for Rational
    end
    @test !(@inferred X ≈ translate(X, N[1 // 1000]))  # above default tolerance for all types
    @test !(@inferred X ≈ S2) && !(@inferred S2 ≈ X) && !(@inferred X ≈ B) &&
          !(@inferred B ≈ X)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(X, S2)
    @test_throws DimensionMismatch isdisjoint(S2, X)
    # disjoint
    Y = Interval(N(3), N(4))
    for (Z, W) in ((X, Y), (Y, X))
        @test @inferred isdisjoint(Z, W)
        @test_broken @inferred isdisjoint(Z, W, true)  # TODO make this type-stable
        res, w = isdisjoint(Z, W, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    # overlapping
    Y = Interval(N(1), N(3))
    for (Z, W) in ((X, X), (X, Y), (Y, X))
        @test !(@inferred isdisjoint(Z, W))
        res, w = isdisjoint(Z, W, true)
        @test !res && w isa Vector{N} && w ∈ Z && w ∈ W
    end
    # tolerance
    if N == Float64
        Y = Interval(2.0 + 1e-9, 3.0)
        @test !(@inferred isdisjoint(X, Y))
        res, w = isdisjoint(X, Y, true)
        # TODO ∈ and isdisjoint should be consistent
        @test_broken !res && w isa Vector{N} && w ∈ X && w ∈ Y
        r = _rtol(N)
        @assert r > N(1e-10) "default tolerance changed; adapt test"
        set_rtol(N, N(1e-10))
        @test @inferred isdisjoint(X, Y)
        res, w = isdisjoint(X, Y, true)
        @test res && w isa Vector{N} && isempty(w)
        # restore tolerance
        set_rtol(N, r)
    end

    # isequal
    @test @inferred X == X
    @test (@inferred X != S2) && (@inferred S2 != X) && (@inferred X != B) && @inferred B != X

    # isequivalent
    @test_throws DimensionMismatch isequivalent(X, S2)
    @test_throws DimensionMismatch isequivalent(S2, X)
    @test @inferred isequivalent(X, X)
    @test !(@inferred isequivalent(X, Interval(N(1), N(2))))
    @test_broken @inferred isequivalent(X, B)  # TODO make this type-stable
    @test isequivalent(X, B) && isequivalent(B, X)

    # isstrictsubset
    @test_throws DimensionMismatch X ⊂ S2
    for Y in (X, B)
        @test_broken @inferred Y ⊂ X  # TODO make this type-stable
        @test !(Y ⊂ X)
        @test_broken @inferred ⊂(Y, X, true)  # TODO make this type-stable
        res, w = ⊂(Y, X, true)
        @test !res && w isa Vector{N} && isempty(w)
    end
    for Y in (Interval(N(-1), N(2)), Interval(N(0), N(3)))
        @test !(Y ⊂ X)
        res, w = ⊂(Y, X, true)
        @test !res && w isa Vector{N} && w ∈ Y && w ∉ X
    end
    for Y in (Interval(N(-1), N(2)), Interval(N(0), N(3)))
        @test_broken @inferred ⊂(X, Y, true)  # TODO make this type-stable
        @test X ⊂ Y
        res, w = ⊂(X, Y, true)
        @test res && w isa Vector{N} && w ∉ X && w ∈ Y
    end

    # issubset
    @test_throws DimensionMismatch X ⊆ S2
    @test_throws DimensionMismatch S2 ⊆ X
    for Y in (X, B)
        @test_broken @inferred X ⊆ Y  # TODO make this type-stable
        @test X ⊆ Y
        @test_broken @inferred ⊆(X, Y, true)  # TODO make this type-stable
        res, w = ⊆(X, Y, true)
        @test res && w isa Vector{N} && w == N[]
    end
    for Y in (Interval(N(0), N(1)), Interval(N(1), N(3)))
        @test X ⊈ Y
        res, w = ⊆(X, Y, true)
        @test !res && w isa Vector{N} && w ∈ X && w ∉ Y
    end

    # linear_combination
    @test_throws DimensionMismatch linear_combination(X, S2)
    @test_throws DimensionMismatch linear_combination(S2, X)
    for Y in (linear_combination(X, Xnc), linear_combination(Xnc, X))
        @test isidentical(Y, Interval(N(0), N(5)))
    end
    Z = @inferred linear_combination(X, X)
    @test isidentical(Z, X)
    @test_broken @inferred linear_combination(X, B)  # TODO make this type-stable
    for Z in (linear_combination(X, B), linear_combination(B, X))
        @test isidentical(Z, X)
    end
    Y = Interval(N(3), N(4))
    for Z in ((@inferred linear_combination(X, Y)), @inferred linear_combination(Y, X))
        @test isidentical(Z, Interval(N(0), N(4)))
    end

    # minkowski_difference
    @test_throws DimensionMismatch minkowski_difference(X, S2)
    @test_throws DimensionMismatch minkowski_difference(S2, X)
    # empty difference
    Y = Interval(N(0), N(3))
    Z = minkowski_difference(X, Y)
    @test Z isa EmptySet{N} && Z == EmptySet{N}(1)
    # nonempty difference
    Y = minkowski_difference(X, X)
    @test isidentical(Y, Interval(N(0), N(0)))
    Y = Interval(N(1), N(3))
    Z = minkowski_difference(X, Y)
    @test isidentical(Z, Interval(N(-1), N(-1)))
    # equivalent sets
    for Y in (minkowski_difference(X, B), minkowski_difference(B, X))
        @test isequivalent(Y, Interval(N(0), N(0)))
    end

    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(X, S2)
    @test_throws DimensionMismatch minkowski_sum(S2, X)
    # Interval + Interval = Interval
    Y = @inferred minkowski_sum(X, X)
    Z = Interval(N(0), N(4))
    @test isidentical(Y, Z)
    # general
    for Y in ((@inferred minkowski_sum(X, B)), @inferred minkowski_sum(B, X))
        @test Y isa LazySet{N} && isequivalent(Y, Z)
    end
end

for N in @tN([Float64, Float32])
    X = Interval(N(0), N(2))

    # rand
    Y = rand(Interval; N=N)
    @test Y isa Interval{N} && dim(Y) == 1
    @test_throws AssertionError rand(Interval; N=N, dim=2)

    # rationalize
    Y = @inferred rationalize(X)
    @test Y isa Interval{Rational{Int64}} && Y == Interval(N(0 // 1), N(2 // 1))

    # exponential_map
    Y = @inferred exponential_map(ones(N, 1, 1), X)
    @test isidentical(Y, Interval(N(0), N(exp(1)) * N(2)))
end
