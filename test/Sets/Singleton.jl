using LazySets, Test
using LazySets.ReachabilityBase.Arrays: ispermutation, SingleEntryVector
using LazySets.ReachabilityBase.Comparison: set_rtol, _rtol
@static if isdefined(Main, :GeometryBasics)
    import GeometryBasics
end
@static if isdefined(Main, :Polyhedra)
    import Polyhedra
end
@static if isdefined(Main, :MiniQhull)
    import MiniQhull
end
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

function isidentical(::Singleton, ::Singleton)
    return false
end

function isidentical(S1::Singleton{N}, S2::Singleton{N}) where {N}
    return S1.element == S2.element
end

for N in @tN([Float64, Float32, Rational{Int}])
    # auxiliary sets
    P = VPolygon([N[2, -1]])  # equivalent set
    Xnc = UnionSet(P, BallInf(N[4, 0], N(1)))  # nonconvex set

    # constructor
    e = N[2, -1]
    S = @inferred Singleton(e)
    S1 = @inferred Singleton(N[2])  # 1D hyperrectangle
    S3 = @inferred Singleton(N[2, -1, 3])  # 3D hyperrectangle
    # constructor from scalars
    S2 = @inferred Singleton(N(2), N(-1))
    @test isidentical(S2, S)
    # constructor from non-Vector
    e2 = @inferred SingleEntryVector(1, 1000, N(2))
    S2 = @inferred Singleton(e2)
    @test element(S2) == e2

    # convert
    # from CartesianProduct of AbstractSingletons
    X = @inferred Singleton(N[2]) × ZeroSet{N}(1)
    @test isidentical(convert(Singleton, X), Singleton(N[2, 0]))

    # an_element
    x = @inferred an_element(S)
    @test x isa Vector{N} && x ∈ S

    # area
    @test_throws DimensionMismatch area(S1)
    res = @inferred area(S)
    @test res === N(0)
    res = @inferred area(S3)
    @test res === N(0)

    # center
    c = @inferred center(S)
    @test c isa AbstractVector{N} && c == e

    # chebyshev_center_radius
    c, r = @inferred chebyshev_center_radius(S)
    @test c isa Vector{N} && c == e && r === N(0)

    # complement
    X = @inferred complement(S)
    @test_broken X isa Universe{N} && dim(X) == 2  # TODO this should maybe change
    @test X isa UnionSetArray{N} && dim(X) == 2

    # concretize
    S2 = @inferred concretize(S)
    @test isidentical(S2, S)

    # constrained_dimensions
    @test (@inferred constrained_dimensions(S)) == 1:2

    # constraints_list
    clist = [HalfSpace(N[1, 0], N(2)), HalfSpace(N[0, 1], N(-1)),
             HalfSpace(N[-1, 0], N(-2)), HalfSpace(N[0, -1], N(1))]
    @test ispermutation((@inferred constraints_list(S)), clist)
    clist2 = constraints_list(S; min_constraints=true)
    @test length(clist2) == 3 && isequivalent(HPolygon(clist), HPolygon(clist2))

    # constraints
    @test ispermutation(collect(@inferred constraints(S)), clist)

    # convex_hull (unary)
    S2 = @inferred convex_hull(S)
    @test isidentical(S, S2)

    # copy
    @test_broken @inferred copy(S)  # TODO make this type-stable
    S2 = copy(S)
    @test isidentical(S, S2)

    # diameter
    @test_throws ArgumentError diameter(S, N(1 // 2))
    for res in ((@inferred diameter(S)), (@inferred diameter(S, Inf)), @inferred diameter(S, 2))
        if N <: AbstractFloat
            @test res === N(0)
        else
            @test res == 0.0
        end
    end

    # dim
    @test @inferred dim(S) == 2
    @test @inferred dim(S3) == 3

    # element
    @test (@inferred element(S)) == e
    for i in 1:2
        @test (@inferred element(S, i)) == e[i]
    end

    # eltype
    @test (@inferred eltype(S)) == N
    @test (@inferred eltype(typeof(S))) == N

    # extrema
    res = @inferred extrema(S)
    @test res isa Tuple{Vector{N},Vector{N}} && res[1] == e && res[2] == e
    @test_throws DimensionMismatch extrema(S, 3)
    res = @inferred extrema(S, 1)
    @test res isa Tuple{N,N} && res[1] == e[1] && res[2] == e[1]

    # generators
    gens = collect(@inferred generators(S))
    @test gens isa AbstractVector{<:AbstractVector{N}} && isempty(gens)

    # genmat
    gens = @inferred genmat(S)
    @test gens isa Matrix{N} && gens == Matrix{N}(undef, 2, 0)

    # high
    res = @inferred high(S)
    @test res isa Vector{N} && res == e
    @test_throws DimensionMismatch high(S, 3)
    res = @inferred high(S, 1)
    @test res === e[1]

    # isbounded
    @test @inferred isbounded(S)

    # isboundedtype
    @test @inferred isboundedtype(typeof(S))

    # isconvex
    @test @inferred isconvex(S)

    # isconvextype
    @test @inferred isconvextype(typeof(S))

    # isempty
    @test !(@inferred isempty(S))
    @test_broken @inferred isempty(S, true)  # TODO make this type-stable (witness)
    res, w = isempty(S, true)
    @test !res && w isa Vector{N} && w ∈ S

    # isflat
    @test @inferred isflat(S)

    # isoperation
    @test !(@inferred isoperation(S))

    # isoperationtype
    @test !(@inferred isoperationtype(typeof(S)))

    # ispolyhedral
    @test @inferred ispolyhedral(S)

    # ispolyhedraltype
    @test @inferred ispolyhedraltype(typeof(S))

    # ispolytopic
    @test @inferred ispolytopic(S)

    # ispolytopictype
    @test @inferred ispolytopictype(typeof(S))

    # isuniversal
    @test !(@inferred isuniversal(S))
    @test_broken @inferred isuniversal(S, true)  # TODO make this type-stable (witness)
    res, w = isuniversal(S, true)
    @test !res && w isa Vector{N} && w ∉ S

    # low
    res = @inferred low(S)
    @test res isa Vector{N} && res == e
    @test_throws DimensionMismatch low(S, 3)
    res = @inferred low(S, 1)
    @test res === e[1]

    # ngens
    @test (@inferred ngens(S)) == 0

    # norm
    @test_throws ArgumentError norm(S, N(1 // 2))
    for res in ((@inferred norm(S)), @inferred norm(S, Inf))
        if N <: AbstractFloat
            @test res === N(2)
        else
            @test res == 2.0
        end
    end
    res = @inferred norm(S, 2)
    if N <: AbstractFloat
        @test res isa N
    end
    @test res == norm(e, 2)

    # radius
    @test_throws ArgumentError radius(S, N(1 // 2))
    for res in ((@inferred radius(S)), (@inferred radius(S, Inf)), @inferred radius(S, 2))
        if N <: AbstractFloat
            @test res === N(0)
        else
            @test res == 0.0
        end
    end

    # radius_hyperrectangle
    res = @inferred radius_hyperrectangle(S)
    @test res isa AbstractVector{N} && res == N[0, 0]
    @test_throws DimensionMismatch radius_hyperrectangle(S, 3)
    res = @inferred radius_hyperrectangle(S, 1)
    @test res === N(0)

    # rectify
    S2 = @inferred rectify(S)
    @test isidentical(S2, Singleton(N[2, 0]))

    # reflect
    S2 = @inferred reflect(S)
    @test isidentical(S2, Singleton(N[-2, 1]))

    # remove_redundant_generators
    @test isidentical((@inferred remove_redundant_generators(S)), S)

    # singleton_list
    res = @inferred singleton_list(S)
    @test res isa Vector{Singleton{N,Vector{N}}}
    @test res == [S]

    # togrep
    Z = @inferred togrep(S)
    @test Z isa Zonotope{N} && isequivalent(Z, S)

    # tosimplehrep
    A, b = @inferred tosimplehrep(S)
    @test A isa Matrix{N} && b isa Vector{N} && size(A) == (4, 2) && length(b) == 4
    # no precise test here; the `tosimplehrep` implementation is tested elsewhere

    # triangulate_faces
    @test_throws DimensionMismatch triangulate_faces(S)
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :GeometryBasics)
        res = @inferred triangulate_faces(S3)
        @test res isa Tuple{Matrix{Float32},Vector{Tuple{Int,Int,Int}}}
        @test length(res[2]) == 0
    end

    # vertices_list
    vlist = @inferred vertices_list(S)
    @test vlist isa Vector{Vector{N}} && vlist == [e]

    # vertices
    res = collect(@inferred vertices(S))
    @test res isa Vector{Vector{N}} && res == vlist

    # volume
    res = @inferred volume(S)
    @test res === N(0)
    res = @inferred volume(S3)
    @test res === N(0)

    # affine_map
    @test_throws DimensionMismatch affine_map(ones(N, 2, 1), S, N[1, 1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 2), S, N[1])
    @test_broken @inferred affine_map(ones(N, 1, 2), S, N[1])  # TODO make this type-stable
    S2 = affine_map(ones(N, 1, 2), S, N[1])
    @test isidentical(S2, Singleton(N[2]))
    @test_broken @inferred affine_map(ones(N, 3, 2), S, N[1, 2, 3])  # TODO make this type-stable
    S2 = affine_map(ones(N, 3, 2), S, N[1, 2, 3])
    @test isidentical(S2, Singleton(N[2, 3, 4]))

    # distance (between point and set)
    @test_throws DimensionMismatch distance(S, N[0])
    @test_throws ArgumentError distance(S, N[0, 0]; p=N(1 // 2))
    x = N[1, 3]
    @test (@inferred distance(x, S)) == (@inferred distance(S, x)) == distance(x, e)
    for p in N[1, 2, Inf]
        @test (@inferred distance(x, S; p=p)) == (@inferred distance(S, x; p=p)) ==
              @inferred distance(x, e; p=p)
    end

    # exponential_map
    @test_throws DimensionMismatch exponential_map(ones(N, 1, 1), S)
    @test_throws DimensionMismatch exponential_map(ones(N, 2, 1), S)
    if N <: AbstractFloat
        S2 = @inferred exponential_map(N[1 0; 0 2], S)
        @test S2 isa Singleton && element(S2) == N[2 * exp(1), -exp(2)]
    end

    # in
    @test_throws DimensionMismatch N[0] ∈ S
    @test (@inferred e ∈ S) && (@inferred N[2, -2] ∉ S)

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0], S)
    @test_throws ArgumentError is_interior_point(N[0, 0], S; ε=N(0))
    @test_throws ArgumentError is_interior_point(N[0, 0], S; p=N(1 // 2))
    if N <: AbstractFloat
        @test !(@inferred is_interior_point(N[2, -1], S))
        @test !(@inferred is_interior_point(N[2, -2], S))
    else
        @test_throws ArgumentError is_interior_point(N[2, -1], S)
        @test !(@inferred is_interior_point(N[2, -1], S; ε=1 // 100))
        @test !(@inferred is_interior_point(N[2, -2], S; ε=1 // 100))
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([1.0, 1], S)
    end

    # linear_map
    @test_throws DimensionMismatch linear_map(ones(N, 1, 1), S)
    S2 = @inferred linear_map(N[1 0; 0 2], S)
    @test isidentical(S2, Singleton(N[2, -2]))
    S2 = @inferred linear_map(zeros(N, 1, 2), S)  # zero map
    @test isidentical(S2, Singleton(N[0]))

    # linear_map_inverse
    @test_throws AssertionError LazySets.linear_map_inverse(ones(N, 1, 1), S)
    # invertible map
    X = LazySets.linear_map_inverse(N[1 0; 0 1], S)
    @test X isa Zonotope{N} && isequivalent(X, S)
    # noninvertible map
    M = ones(N, 2, 1)
    X = LazySets.linear_map_inverse(M, S)
    A, b = tosimplehrep(S)
    Y = HPolytope(A * M, b)
    @test X isa LazySet{N}
    @test_broken isequivalent(X, Y)  # TODO this should work

    # permute
    @test_throws DimensionMismatch permute(S, [1])
    @test_throws DimensionMismatch permute(S, [1, -1])
    @test_throws DimensionMismatch permute(S, [1, 3])
    @test_throws ArgumentError permute(S, [1, 1])
    S2 = @inferred permute(S, [1, 2])
    @test isidentical(S2, S)
    S2 = @inferred permute(S, [2, 1])
    @test isidentical(S2, Singleton(N[-1, 2]))

    # project
    @test_throws DimensionMismatch project(S, [1, 2, 3])
    @test_throws DimensionMismatch project(S, [1, -1])
    @test_throws DimensionMismatch project(S, [1, 3])
    @test_throws ArgumentError project(S, [1, 1])
    S2 = @inferred project(S, [1])
    @test isidentical(S2, Singleton(N[2]))
    S2 = @inferred project(Singleton(N[4, 3, 2, 1]), [2, 4])
    @test isidentical(S2, Singleton(N[3, 1]))

    # sample
    res = @inferred sample(S)
    @test res isa Vector{N} && res ∈ S
    res = @inferred sample(S, 2)
    @test res isa Vector{Vector{N}} && length(res) == 2 && all(x ∈ S for x in res)

    # scale
    S2 = @inferred scale(N(2), S)
    @test isidentical(S2, Singleton(N[4, -2]))
    # degenerate case
    S2 = @inferred scale(N(0), S)
    @test isidentical(S2, Singleton(N[0, 0]))
    # scale!
    S2 = copy(S)
    @inferred scale!(N(2), S2)
    @test isidentical(S2, Singleton(N[4, -2]))
    # degenerate case
    S2 = copy(S)
    @inferred scale!(N(0), S2)
    @test isidentical(S2, Singleton(N[0, 0]))
    # negative factor
    S2 = copy(S)
    @inferred scale!(N(-2), S2)
    @test isidentical((@inferred scale(N(-2), S)), S2) && isidentical(S2, Singleton(N[-4, 2]))

    # split
    @test_throws ArgumentError split(S, [0, 2])
    X = Hyperrectangle(N[2, -1], N[0, 0])
    @test (@inferred split(S, [2, 2])) == [X, X, X, X]
    # split at given bounds
    @test (@inferred split(S, [N[2], N[-1]])) == [X, X, X, X]

    # support_function
    @test_throws DimensionMismatch ρ(N[1], S)
    res = @inferred ρ(N[1, 1], S)
    @test res === N(1)
    res = @inferred ρ(N[-1, -1], S)
    @test res === N(-1)

    # support_vector
    @test_throws DimensionMismatch σ(N[1], S)
    res = @inferred σ(N[1, 1], S)
    @test res isa Vector{N} && res == e
    res = @inferred σ(N[-1, -1], S)
    @test res isa Vector{N} && res == e

    # translate
    @test_throws DimensionMismatch translate(S, N[1])
    @test_broken @inferred translate(S, N[1, 2])  # TODO make this type-stable
    S2 = translate(S, N[1, 2])
    @test isidentical(S2, Singleton(N[3, 1]))
    # translate!
    @test_throws DimensionMismatch translate!(S, N[1])
    S2 = copy(S)
    @inferred translate!(S2, N[1, 2])
    @test isidentical(S2, Singleton(N[3, 1]))

    # cartesian_product
    S2a = Singleton(N[2])
    S2b = Singleton(N[-1])
    S2 = @inferred cartesian_product(S2a, S2b)
    @test isidentical(S2, S)

    # convex_hull (binary)
    @test_throws DimensionMismatch convex_hull(S, S3)
    X = convex_hull(S, S)
    @test X isa LazySet{N} && isequivalent(X, S)
    e2 = N[3, 1]
    S2 = Singleton(e2)
    Y = LineSegment(e, e2)
    for X in (convex_hull(S, S2), convex_hull(S2, S))
        @test X isa LazySet{N} && isequivalent(X, Y)
    end

    # difference
    @test_throws DimensionMismatch difference(S, S3)
    @static if isdefined(@__MODULE__, :IntervalBoxes)
        # full set
        X = difference(S, Singleton(N[1, 2]))
        @test isequivalent(X, S)
        # empty set
        X = difference(S, S)
        @test isequivalent(X, EmptySet{N}(2))
    end

    # distance (between two sets)
    @test_throws DimensionMismatch distance(S, S3)
    @test_throws ArgumentError distance(S, S; p=N(1 // 2))
    res = @inferred distance(S, S)
    if N <: AbstractFloat
        @test res === N(0)
    else
        @test res isa Float64 && res == 0.0
    end
    e2 = N[3, 1]
    S2 = Singleton(e2)
    @test (@inferred distance(S, S2)) == (@inferred distance(S2, S)) == @inferred distance(e, e2)
    for p in N[1, 2, Inf]
        @test (@inferred distance(S, S2; p=p)) == (@inferred distance(S2, S; p=p)) ==
              @inferred distance(e, e2; p=p)
    end

    # exact_sum
    @test_throws DimensionMismatch exact_sum(S, S3)
    S2 = Singleton(N[3, 1])
    for S2b in ((@inferred exact_sum(S, S2)), @inferred exact_sum(S2, S))
        @test isidentical(S2b, Singleton(N[5, 0]))
    end

    # intersection
    @test_throws DimensionMismatch intersection(S, S3)
    # disjoint
    X = intersection(S, Singleton(N[3, 1]))
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    # overlapping
    S2 = intersection(S, S)
    @test isidentical(S2, S)

    # isapprox
    @test @inferred S ≈ S
    res = @inferred S ≈ translate(S, N[1 // 100000000, 0])
    if N <: AbstractFloat
        @test res  # below default tolerance for AbstractFloat
    else
        @test !res  # zero default tolerance for Rational
    end
    @test !(@inferred S ≈ translate(S, N[1 // 1000, 0]))  # above default tolerance for all types
    @test !(@inferred S ≈ S3) && !(@inferred S3 ≈ S) && !(@inferred S ≈ P) && !(@inferred P ≈ S)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(S, S3)
    # disjoint
    S2 = Singleton(N[3, 1])
    for (S2a, S2b) in ((S, S2), (S2, S))
        if (N == Float64 && VERSION < v"1.12") || (N == Rational{Int} && VERSION < v"1.11")
            @test @inferred isdisjoint(S2a, S2b)
        else
            @test_broken @inferred isdisjoint(S2a, S2b)  # TODO make this type-stable (witness)
            @test isdisjoint(S2a, S2b)
        end
        @test_broken @inferred isdisjoint(S2a, S2b, true)  # TODO make this type-stable (witness)
        res, w = isdisjoint(S2a, S2b, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    # overlapping
    @test !isdisjoint(S, S)
    res, w = isdisjoint(S, S, true)
    @test !res && w isa Vector{N} && w ∈ S
    # tolerance
    if N == Float64
        S2 = Singleton(N[2 + 1e-9, -1])
        @test !isdisjoint(S, S2)
        r = _rtol(N)
        @assert r > N(1e-10) "default tolerance changed; adapt test"
        set_rtol(N, N(1e-10))
        @test isdisjoint(S, S2)
        # restore tolerance
        set_rtol(N, r)
    end

    # isequal
    @test @inferred S == S
    @test (@inferred S != S3) && (@inferred S3 != S) && (@inferred S != P) && @inferred P != S

    # isequivalent
    @test_throws DimensionMismatch isequivalent(S, S3)
    @test @inferred isequivalent(S, S)
    @test !(@inferred isequivalent(S, Singleton(N[3, 1])))
    @test_broken (@inferred isequivalent(S, P)) && @inferred isequivalent(P, S)  # TODO make this type-stable (witness)
    @test isequivalent(S, P) && isequivalent(P, S)

    # isstrictsubset
    @test_throws DimensionMismatch S ⊂ S3
    for X in (S, P)
        @test_broken !(@inferred X ⊂ S)  # TODO make this type-stable (witness)
        @test !(X ⊂ S)
        res, w = ⊂(X, S, true)
        @test !res && w isa Vector{N} && isempty(w)
    end
    X = EmptySet{N}(2)
    @test X ⊂ S
    res, w = ⊂(X, S, true)
    @test res && w isa Vector{N} && w ∉ X && w ∈ S

    # issubset
    @test_throws DimensionMismatch S ⊆ S3
    if (N == Float64 && VERSION < v"1.12") || (N == Rational{Int} && VERSION < v"1.11")
        @test @inferred S ⊆ S
    else
        @test_broken @inferred S ⊆ S  # TODO make this type-stable (witness)
        @test S ⊆ S
    end
    @test_broken @inferred ⊆(S, S, true)  # TODO make this type-stable (witness)
    res, w = ⊆(S, S, true)
    @test res && w isa Vector{N} && w == N[]
    if VERSION >= v"1.12"
        @test_broken @inferred S ⊆ P  # TODO make this type-stable (witness)
        @test S ⊆ P
    else
        @test @inferred S ⊆ P
    end
    res, w = ⊆(S, P, true)
    @test res && w isa Vector{N} && w == N[]
    S2 = Singleton(N[3, 1])
    @test S ⊈ S2
    res, w = ⊆(S, S2, true)
    @test !res && w isa Vector{N} && w ∈ S && w ∉ S2

    # linear_combination
    @test_throws DimensionMismatch linear_combination(S, S3)
    @test_broken linear_combination(S, Xnc) isa LazySet{N}  # TODO implement `linear_combination` for non-convex sets
    @test_broken linear_combination(Xnc, S) isa LazySet{N}
    S2 = linear_combination(S, S)
    @test isidentical(S2, S)
    for X in (linear_combination(S, P), linear_combination(P, S))
        @test X isa LazySet{N} && isequivalent(X, S)
    end

    # minkowski_difference
    @test_throws DimensionMismatch minkowski_difference(S, S3)
    # empty difference
    X = minkowski_difference(S, Hyperrectangle(N[2, -1], N[1, 1]))
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    # equivalent sets: only the origin remains
    @test_broken @inferred minkowski_difference(S, S)  # TODO make this type-stable
    S2 = minkowski_difference(S, S)
    @test isidentical(S2, Singleton(N[0, 0]))
    for X in ((@inferred minkowski_difference(S, P)), @inferred minkowski_difference(P, S))
        vlist = vertices_list(X)
        @test length(vlist) == 1 && iszero(vlist)
    end
    # nonempty difference
    @test_broken @inferred minkowski_difference(S, Singleton(N[3, 1]))  # TODO make this type-stable
    S2 = minkowski_difference(S, Singleton(N[3, 1]))
    @test isidentical(S2, Singleton(N[-1, -2]))

    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(S, S3)
    S2 = Singleton(N[3, 1])
    for S2b in ((@inferred minkowski_sum(S, S2)), @inferred minkowski_sum(S2, S))
        @test isidentical(S2b, Singleton(N[5, 0]))
    end
end

for N in @tN([Float64, Float32])
    S = Singleton(N[2, -1])

    # rand
    S2 = rand(Singleton; N=N)
    @test S2 isa Singleton{N} && dim(S2) == 2
    S2 = rand(Singleton; N=N, dim=3)
    @test S2 isa Singleton{N} && dim(S2) == 3

    # rationalize
    S2 = rationalize(S)
    T = Rational{Int64}
    @test S2 isa Singleton{T} && isidentical(S2, Singleton(T[2, -1]))
end

for N in @tN([Float64, Rational{Int}])
    S = Singleton(N[2, -1])

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        X = polyhedron(S)
        @test X isa Polyhedra.DefaultPolyhedron{N}
        @test ispermutation(constraints_list(convert(HPolytope, X)), constraints_list(S))
    end

    # triangulate
    @static if isdefined(@__MODULE__, :MiniQhull)
        @test_broken triangulate(S)  # TODO this should work
    end
end
