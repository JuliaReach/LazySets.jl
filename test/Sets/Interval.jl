using LazySets, Test
import IntervalArithmetic as IA
@static if VERSION >= v"1.9"
    vIA = pkgversion(IA)
else
    import PkgVersion
    vIA = PkgVersion.Version(IA)
end
using LazySets.ReachabilityBase.Arrays: ispermutation
using LazySets.ReachabilityBase.Arrays: SingleEntryVector

function isidentical(::Interval, ::Interval)
    return false
end

function isidentical(X1::Interval{N}, X2::Interval{N}) where {N}
    return X1.dat == X2.dat
end

for N in [Float64, Float32, Rational{Int}]
    # auxiliary sets
    X2 = Singleton(N[0, 0])  # 2D set
    B = BallInf(N[1], N(1))  # equivalent set
    PZ = SimpleSparsePolynomialZonotope(N[0], zeros(N, 1, 0), zeros(Int, 0, 0))  # nonconvex

    # default constructor from IntervalArithmetic.Interval
    itv = IA.interval(N(0), N(2))
    X = Interval(itv)
    @test X isa Interval{N} && X.dat == itv

    # constructors from two numbers, from a vector, and with promotion
    for Y in (Interval(N(0), N(2)), Interval(N[0, 2]), Interval(0, N(2)))
        @test isidentical(Y, X)
    end

    # constructor from a number
    X0 = Interval(N(0))
    @test X0 isa Interval{N} && X0 == Interval(N(0), N(0))

    # unbounded intervals are caught
    @test_throws AssertionError Interval(-Inf, zero(N))
    @test_throws AssertionError Interval(zero(N), Inf)

    # convert
    # to and from IntervalArithmetic.Interval
    Y = convert(IA.Interval, X)
    @test Y == X.dat
    Z = convert(Interval, Y)
    @test isidentical(Z, X)
    # from hyperrectangular set
    @test_throws AssertionError convert(Interval, Hyperrectangle(N[0, 0], N[1, 1]))
    Y = Hyperrectangle(N[1], N[1])
    Z = convert(Interval, Y)
    @test isidentical(Z, X)
    # from Rectification of Interval
    Y = Interval(N(-1), N(2))
    Z = convert(Interval, Rectification(Y))
    @test isidentical(Z, X)
    # from MinkowskiSum of Intervals
    Z = convert(Interval, Y + Y)
    @test isidentical(Z, Interval(N(-2), N(4)))
    # from LazySet
    M = hcat(N[1])
    Y = convert(Interval, M * X)
    @test isidentical(Y, X)
    # from non-convex set
    Y = UnionSet(Interval(1, 2), Interval(3, 4))
    @test_throws AssertionError convert(Interval, Y)

    # an_element
    x = an_element(X)
    @test x isa Vector{N} && length(x) == 1 && x[1] isa N && N(0) <= x[1] <= N(2)

    # area
    @test_throws DimensionMismatch area(X)

    # center
    c = center(X)
    @test c isa Vector{N} && c == [N(1)]
    v = center(X, 1)
    @test v isa N && v == N(1)
    @test_throws AssertionError center(X, 2)

    # chebyshev_center_radius
    c, r = chebyshev_center_radius(X)
    @test c isa Vector{N} && c == [N(1)]
    @test r isa N && r == N(1)

    # complement
    Y = complement(X)
    @test Y isa UnionSet{N}
    @test ispermutation(array(Y), [HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-2))])

    # concretize
    Y = concretize(X)
    @test isidentical(Y, X)

    # constrained_dimensions
    @test constrained_dimensions(X) == 1:1

    # constraints_list
    cs = constraints_list(X)
    @test ispermutation(cs, [HalfSpace(N[1], N(2)), HalfSpace(N[-1], N(0))])

    # constraints
    cs2 = collect(constraints(X))
    @test ispermutation(cs2, cs)

    # convex_hull (unary)
    Y = convex_hull(X)
    @test isidentical(Y, X)

    # copy
    Y = copy(X)
    @test isidentical(Y, X)

    # diameter
    for res in (diameter(X), diameter(X, Inf), diameter(X, 2))
        @test res isa N && res == N(2)
    end

    # dim
    @test dim(X) == 1

    # eltype
    @test eltype(X) == N
    @test eltype(typeof(X)) == N

    # extrema
    res = extrema(X)
    @test res isa Tuple{Vector{N},Vector{N}} && res[1] == N[0] && res[2] == N[2]
    res = extrema(X, 1)
    @test res isa Tuple{N,N} && res[1] == N(0) && res[2] == N(2)
    @test_throws AssertionError extrema(X, 2)

    # generators
    @test collect(generators(X)) == [N[1]]
    # degenerate case
    gens = collect(generators(X0))
    @test gens isa Vector{SingleEntryVector{N}} && isempty(gens)

    # genmat
    @test genmat(X) == hcat(N[1])
    # degenerate case
    gens = genmat(X0)
    @test gens isa Matrix{N} && gens == Matrix{N}(undef, 1, 0)

    # high
    res = high(X)
    @test res isa Vector{N} && res == N[2]
    res = high(X, 1)
    @test res isa N && res == N(2)
    @test_throws AssertionError high(X, 2)

    # isbounded
    @test isbounded(X)

    # isboundedtype
    @test isboundedtype(typeof(X))

    # isconvextype
    @test isconvextype(typeof(X))

    # isempty
    @test !isempty(X)
    res, w = isempty(X, true)
    @test !res && w ∈ X
    if N == Float64
        @test w isa Vector{N}
    else
        @test_broken w isa Vector{N}  # TODO this should change
    end

    # isflat
    ztol = LazySets._ztol(N)  # default absolute zero tolerance
    @test isflat(Interval(N(0), ztol))
    @test !isflat(Interval(N(0), 2 * ztol + N(1 // 100)))

    # isoperation
    @test !isoperation(X)

    # isoperationtype
    @test !isoperationtype(typeof(X))

    # ispolyhedral
    @test ispolyhedral(X)

    # isuniversal
    @test !isuniversal(X)
    res, w = isuniversal(X, true)
    @test !res && w isa Vector{N} && w ∉ X

    # low
    res = low(X)
    @test res isa Vector{N} && res == N[0]
    res = low(X, 1)
    @test res isa N && res == N(0)
    @test_throws AssertionError low(X, 2)

    # min
    v = min(X)
    @test v isa N && v == N(0)

    # max
    v = max(X)
    @test v isa N && v == N(2)

    # ngens
    @test ngens(X) == 1
    # degenerate case
    @test ngens(X0) == 0

    # norm
    for Y in (X, Interval(N(-2), N(1)))
        for res in (norm(Y), norm(Y, Inf), norm(Y, 2))
            @test res isa N && res == N(2)
        end
    end

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        Y = polyhedron(X)
        @test Y isa Polyhedra.Interval{N} && ispermutation(Y.vrep.points.points, [[N(0)], [N(2)]])
    end

    # radius
    for res in (radius(X), radius(X, Inf), radius(X, 2))
        @test res isa N && res == N(1)
    end

    # radius_hyperrectangle
    r = radius_hyperrectangle(X)
    @test r isa Vector{N} && r == [N(1)]
    v = radius_hyperrectangle(X, 1)
    @test v isa N && v == N(1)

    # rectify
    Y = Interval(N(-1), N(2))
    Z = rectify(Y)
    @test isidentical(Z, X)
    Y = Interval(N(-2), N(-1))
    Z = rectify(Y)
    @test isidentical(Z, X0)
    Y = Interval(N(1), N(2))
    Z = rectify(Y)
    @test isidentical(Y, Z)

    # reflect
    Y = reflect(X)
    @test isidentical(Y, Interval(N(-2), N(0)))

    # singleton_list
    res = singleton_list(X)
    @test res isa Vector{Singleton{N,Vector{N}}}
    @test ispermutation(res, [Singleton(N[0]), Singleton(N[2])])

    # tosimplehrep
    A, b = tosimplehrep(X)
    @test A isa Matrix{N} && b isa Vector{N}
    @test (A == hcat(N[1, -1]) && b == N[2, 0]) || (A == hcat(N[-1, 1]) && b == N[0, 2])

    # triangulate
    @static if isdefined(@__MODULE__, :MiniQhull)
        Y = triangulate(X)  # TODO does it make sense to triangulate a set with two vertices?
        @test Y isa UnionSetArray{N} && length(Y) == 1
        Y = Y[1]
        @test Y isa VPolytope{N}
        res = vertices_list(Y)
        @test ispermutation(res, [N[0], N[2]])
    end

    # triangulate_faces
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test_throws ArgumentError triangulate_faces(X)
    end

    # vertices_list
    res = vertices_list(X)
    @test res isa Vector{Vector{N}} && ispermutation(res, [N[0], N[2]])
    # degenerate case
    res = vertices_list(X0)
    @test res isa Vector{Vector{N}} && res == [[N(0)]]

    # vertices
    res = collect(vertices(X))
    @test res isa Vector{Vector{N}} && ispermutation(res, [N[0], N[2]])

    # volume
    @test volume(X) == N(2)

    # affine_map
    @test_throws DimensionMismatch affine_map(ones(N, 1, 2), X, N[1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 1), X, N[1])
    Y = affine_map(ones(N, 1, 1), X, N[1])
    @test isidentical(Y, Interval(N(1), N(3)))
    Y = affine_map(ones(N, 2, 1), X, N[1, 2])
    @test Y isa LazySet{N} && isequivalent(Y, LineSegment(N[1, 2], N[3, 4]))

    # distance (between point and set)
    @test_throws DimensionMismatch distance(X, N[0, 0])
    for (x, v) in ((N[1], N(0)), (N[4], N(2)))
        for res in (distance(X, x), distance(x, X))
            @test res == v
            if N <: AbstractFloat
                @test res isa N
            end
        end
    end

    # exponential_map
    @test_throws DimensionMismatch exponential_map(ones(N, 2, 2), X)
    @test_throws DimensionMismatch exponential_map(ones(N, 2, 1), X)

    # in
    @test_throws DimensionMismatch N[0, 0] ∈ X
    @test N[1] ∈ X && N[2] ∈ X && N[3] ∉ X
    # number in interval is invalid
    @test_throws MethodError N(1) ∈ X

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0, 0], X)
    if N <: AbstractFloat
        @test is_interior_point(N[1], X)
        if N == Float64
            @test_broken !is_interior_point(N[2], X)  # TODO fix
        else
            @test !is_interior_point(N[2], X)
        end
        @test !is_interior_point(N[3], X)
    else
        @test_throws ArgumentError is_interior_point(N[1], X)
        @test is_interior_point(N[1], X; ε=1 // 100)
        @test !is_interior_point(N[2], X; ε=1 // 100)
        @test !is_interior_point(N[3], X; ε=1 // 100)
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([1.0], X)
    end

    # linear_map
    @test_throws AssertionError linear_map(ones(N, 2, 2), X)
    Y = linear_map(2 * ones(N, 1, 1), X)
    @test Y isa Interval{N} && isequivalent(Y, Interval(N(0), N(4)))
    Y = linear_map(zeros(N, 1, 1), X)
    @test Y isa Interval{N}
    if vIA == v"0.21.0"
        @test_broken isequivalent(Y, Interval(N(0), N(0)))  # bug in IntervalArithmetic: 0 * I == I
    else
        @test isequivalent(Y, Interval(N(0), N(0)))
    end
    Y = linear_map(ones(N, 2, 1), X)
    @test Y isa LazySet{N} && isequivalent(Y, LineSegment(N[0, 0], N[2, 2]))
    Y = linear_map(zeros(N, 2, 1), X)
    @test Y isa LazySet{N} && isequivalent(Y, ZeroSet{N}(2))

    # linear_map_inverse
    @test_throws AssertionError LazySets.linear_map_inverse(ones(N, 2, 2), X)
    Y = LazySets.linear_map_inverse(ones(N, 1, 1), X)
    @test Y isa LazySet{N} && isequivalent(Y, X)
    Y = LazySets.linear_map_inverse(ones(N, 1, 2), X)
    Z = HPolyhedron([HalfSpace(N[1, 1], N(2)), HalfSpace(N[-1, -1], N(0))])
    @test Y isa LazySet{N} && isequivalent(Y, Z)

    # permute
    @test_throws AssertionError permute(X, [-1])
    @test_throws AssertionError permute(X, [1, 2])
    @test_throws AssertionError permute(X, [2])
    Y = permute(X, [1])
    @test isidentical(Y, X)

    # project
    @test_throws AssertionError project(X, [-1])
    @test_throws AssertionError project(X, [1, 2])
    @test_throws AssertionError project(X, [2])
    Y = project(X, [1])
    @test isidentical(Y, X)

    # sample
    res = sample(X)
    @test res isa Vector{N} && res in X
    res = sample(X, 2)
    @test res isa Vector{Vector{N}} && length(res) == 2 && all(x in X for x in res)

    # scale
    Y = scale(N(2), X)
    @test isidentical(Y, Interval(N(0), N(4)))
    # degenerate case
    Y = scale(N(0), X)
    @test isidentical(Y, X0)
    # scale!
    @test_throws MethodError scale!(N(2), X)

    # split
    @test_throws AssertionError split(X, 0)
    @test_throws AssertionError split(X, [4, 4])
    Xs = [Interval(N(0), N(1 // 2)), Interval(N(1 // 2), N(1)),
          Interval(N(1), N(3 // 2)), Interval(N(3 // 2), N(2))]
    Ys = split(X, 4)
    @test Ys isa Vector{Interval{N}} && Ys == split(X, [4]) == Xs

    # support_function
    @test_throws AssertionError ρ(N[1, 1], X)
    res = ρ(N[2], X)
    @test res isa N && res == N(4)
    res = ρ(N[-2], X)
    @test res isa N && res == N(0)

    # support_vector
    @test_throws AssertionError σ(N[1, 1], X)
    res = σ(N[2], X)
    @test res isa Vector{N} && res == [N(2)]
    res = σ(N[-2], X)
    @test res isa Vector{N} && res == [N(0)]

    # translate
    @test_throws AssertionError translate(X, N[1, 1])
    Y = translate(X, N[1])
    @test isidentical(Y, Interval(N(1), N(3)))
    # translate!
    @test_throws MethodError translate!(X, N[1, 1])  # TODO this should maybe change
    @test_throws MethodError translate!(X, N[1])  # TODO this should maybe change

    # cartesian_product
    Y = Interval(N(-3), N(1))
    Z = cartesian_product(X, Y)
    @test Z isa Hyperrectangle{N} && Z == Hyperrectangle(N[1, -1], N[1, 2])
    Z = cartesian_product(Y, X)
    @test Z isa Hyperrectangle{N} && Z == Hyperrectangle(N[-1, 1], N[2, 1])

    # convex_hull (binary)
    @test_throws DimensionMismatch convex_hull(X, X2)
    Y = convex_hull(X, X)
    @test isidentical(Y, X)
    Y = Interval(N(-3), N(-1))
    Z = Interval(N(-3), N(2))
    for W in (convex_hull(X, Y), convex_hull(Y, X))
        @test W isa LazySet{N} && isidentical(W, Z)
    end

    # difference
    @test_throws DimensionMismatch difference(X, X2)
    @test_throws DimensionMismatch difference(X2, X)
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
    @test_throws DimensionMismatch distance(X, X2)
    @test_throws DimensionMismatch distance(X2, X)
    for (Y, v) in ((Interval(N(-1), N(1)), N(0)), (Interval(N(4), N(5)), N(2)))
        for res in (distance(X, Y), distance(Y, X))
            @test res isa N && res == v
        end
    end

    # exact_sum
    @test_throws DimensionMismatch exact_sum(X, X2)
    @test_throws DimensionMismatch exact_sum(X2, X)
    Y = Interval(N(3), N(4))
    for Z in (exact_sum(X, Y), exact_sum(Y, X))
        @test isidentical(Z, Interval(N(3), N(6)))
    end

    # intersection
    @test_throws DimensionMismatch intersection(X, X2)
    # disjoint
    Y = intersection(X, Interval(N(3), N(4)))
    @test Y isa EmptySet{N} && Y == EmptySet{N}(1)
    # overlapping
    Y = intersection(X, Interval(N(1), N(3)))
    @test isidentical(Y, Interval(N(1), N(2)))

    # isapprox
    @test X ≈ X
    res = (X ≈ translate(X, N[1 // 100000000]))
    if N <: AbstractFloat
        @test res  # below default tolerance for AbstractFloat
    else
        @test !res  # zero default tolerance for Rational
    end
    @test !(X ≈ translate(X, N[1 // 1000]))  # above default tolerance for all types
    @test !(X ≈ X2) && !(X2 ≈ X) && !(X ≈ B) && !(B ≈ X)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(X, X2)
    # disjoint
    Y = Interval(N(3), N(4))
    @test isdisjoint(X, Y) && isdisjoint(Y, X)
    for (pair, Z) in ((isdisjoint(X, Y, true), Y), (isdisjoint(Y, X, true), Y))
        res, w = pair
        @test res && w isa Vector{N} && isempty(w)
    end
    # overlapping
    Y = Interval(N(1), N(3))
    @test !isdisjoint(X, X) && !isdisjoint(X, Y) && !isdisjoint(Y, X)
    for (pair, Z) in ((isdisjoint(X, X, true), X), (isdisjoint(X, Y, true), Y),
                      (isdisjoint(Y, X, true), Y))
        res, w = pair
        @test !res && w isa Vector{N} && w ∈ X && w ∈ Z
    end
    # tolerance
    if N == Float64
        Y = Interval(2.0 + 1e-9, 3.0)
        @test !isdisjoint(X, Y)
        LazySets.set_rtol(Float64, 1e-10)
        @test isdisjoint(X, Y)
        # restore tolerance
        LazySets.set_rtol(Float64, LazySets.default_tolerance(Float64).rtol)
    end

    # isequal
    @test X == X
    @test X != X2 && X2 != X && X != B && B != X

    # isequivalent
    @test_throws DimensionMismatch isequivalent(X, X2)
    @test_throws DimensionMismatch isequivalent(X2, X)
    @test isequivalent(X, X)
    @test !isequivalent(X, Interval(N(1), N(2)))
    @test isequivalent(X, B) && isequivalent(B, X)

    # isstrictsubset
    @test_throws DimensionMismatch X ⊂ X2
    @test_throws DimensionMismatch X2 ⊂ X
    for Y in (X, B, Interval(N(-1), N(2)), Interval(N(0), N(3)))
        @test !(Y ⊂ X)
        res, w = ⊂(Y, X, true)
        @test !res && w isa Vector{N} && isempty(w)
    end
    for Y in (Interval(N(-1), N(2)), Interval(N(0), N(3)))
        @test X ⊂ Y
        res, w = ⊂(X, Y, true)
        @test res && w isa Vector{N} && w ∉ X && w ∈ Y
    end

    # issubset
    @test_throws DimensionMismatch X ⊆ X2
    @test_throws DimensionMismatch X2 ⊆ X
    for Y in (X, B)
        @test X ⊆ Y
        res, w = ⊆(X, Y, true)
        @test res && w isa Vector{N} && w == N[]
    end
    for Y in (Interval(N(0), N(1)), Interval(N(1), N(3)))
        @test X ⊈ Y
        res, w = ⊆(X, Y, true)
        @test !res && w isa Vector{N} && w ∈ X && w ∉ Y
    end

    # linear_combination
    @test_throws DimensionMismatch linear_combination(X, X2)
    @test_throws ArgumentError linear_combination(X, PZ)
    @test_throws ArgumentError linear_combination(PZ, X)
    for Z in (linear_combination(X, X), linear_combination(X, B), linear_combination(B, X))
        @test isidentical(Z, X)
    end
    Y = Interval(N(3), N(4))
    for Z in (linear_combination(X, Y), linear_combination(Y, X))
        @test isidentical(Z, Interval(N(0), N(4)))
    end

    # minkowski_difference
    @test_throws DimensionMismatch minkowski_difference(X, X2)
    @test_throws DimensionMismatch minkowski_difference(X2, X)
    # empty difference
    Y = Interval(N(0), N(3))
    Z = minkowski_difference(X, Y)
    @test Z isa EmptySet{N} && Z == EmptySet{N}(1)
    # nonempty difference
    Y = minkowski_difference(X, X)
    isidentical(Y, Interval(N(0), N(0)))
    Y = Interval(N(1), N(3))
    Z = minkowski_difference(X, Y)
    @test isidentical(Z, Interval(N(-1), N(-1)))
    for Y in (minkowski_difference(X, B), minkowski_difference(B, X))
        @test isequivalent(Y, Interval(N(0), N(0)))
    end
    # Universe
    U = Universe{N}(1)
    U2 = minkowski_difference(U, X)
    @test U2 isa Universe{N} && dim(U2) == 1

    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(X, X2)
    @test_throws DimensionMismatch minkowski_sum(X2, X)
    # Interval + Interval = Interval
    Y = minkowski_sum(X, X)
    Z = Interval(N(0), N(4))
    @test isidentical(Y, Z)
    # general
    for Y in (minkowski_sum(X, B), minkowski_sum(B, X))
        @test Y isa LazySet{N} && isequivalent(Y, Z)
    end
end

for N in [Float64, Float32]
    X = Interval(N(0), N(2))

    # rand
    Y = rand(Interval; N=N)
    @test Y isa Interval{N} && dim(Y) == 1
    @test_throws AssertionError rand(Interval; N=N, dim=2)

    # rationalize
    Y = rationalize(X)
    @test Y isa Interval{Rational{Int64}} && Y == Interval(N(0 // 1), N(2 // 1))

    # exponential_map
    Y = exponential_map(ones(N, 1, 1), X)
    @test isidentical(Y, Interval(N(0), N(exp(1)) * N(2)))
end
