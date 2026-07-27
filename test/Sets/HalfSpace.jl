using LazySets, Test, SparseArrays
@static if isdefined(Main, :CDDLib)
    import CDDLib
end
@static if isdefined(Main, :Polyhedra)
    import Polyhedra
end
@static if isdefined(Main, :MiniQhull)
    import MiniQhull
end
@static if isdefined(Main, :Symbolics)
    import Symbolics
    using Symbolics: @variables
end
@static if isdefined(Main, :SymEngine)
    import SymEngine
    LazySetsSymEngineExt = Base.get_extension(LazySets, :LazySetsSymEngineExt)
end
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

function isidentical(::HalfSpace, ::HalfSpace)
    return false
end

function isidentical(H1::HalfSpace{N}, H2::HalfSpace{N}) where {N}
    return H1.a == H2.a && H1.b == H2.b
end

for N in @tN([Float64, Float32, Rational{Int}])
    # auxiliary sets
    P = HPolyhedron([HalfSpace(N[1, 2], N(3))])  # equivalent set
    Xnc = UnionSet(P, BallInf(N[3, 0], N(1)))  # nonconvex set

    # constructor
    @test_throws AssertionError HalfSpace(N[0, 0], N(1))  # zero normal vector
    H = @inferred HalfSpace(N[1, 2], N(3))
    @test H isa HalfSpace{N}
    @test H.a == N[1, 2]
    @test H.b == N(3)
    H1 = HalfSpace(N[3], N(4))
    @test H1 isa HalfSpace{N}
    @test H1.a == N[3]
    @test H1.b == N(4)
    H3 = HalfSpace(N[1, 2, 3], N(6))
    @test H3 isa HalfSpace{N}
    @test H3.a == N[1, 2, 3]
    @test H3.b == N(6)

    # convert
    # TODO
    # of normal vector
    H2a = HalfSpace(sparsevec([2], N[1], 3), N(1))
    H2b = convert(HalfSpace{N,Vector{N}}, H2a)
    @test isidentical(H2b, HalfSpace(N[0, 1, 0], N(1))) && H2b.a isa Vector{N}

    # an_element
    x = @inferred an_element(H)
    @test x isa Vector{N} && x ∈ H

    # area
    @test_throws DimensionMismatch area(H1)
    @test_throws AssertionError area(H)  # TODO this should be an ArgumentError

    # chebyshev_center_radius
    @static if isdefined(@__MODULE__, :Polyhedra)  # TODO Polyhedra or some other dep?
        @test_throws ArgumentError chebyshev_center_radius(H)
    end

    # complement
    H2 = @inferred complement(H)
    @test isidentical(H2, HalfSpace(N[-1, -2], N(-3)))

    # concretize
    H2 = @inferred concretize(H)
    @test isidentical(H2, H)

    # constrained_dimensions
    @test (@inferred constrained_dimensions(H)) == 1:2
    @test (@inferred constrained_dimensions(HalfSpace(N[0, 1, 0], N(1)))) == [2]

    # constraints_list
    clist = @inferred constraints_list(H)
    @test clist == [H]

    # constraints
    @test collect(@inferred constraints(H)) == clist

    # convex_hull (unary)
    H2 = @inferred convex_hull(H)
    @test isidentical(H2, H)

    # copy
    @test_broken @inferred copy(H)  # TODO specialize
    H2 = copy(H)
    @test isidentical(H2, H) && H !== H2 && H.a !== H2.a

    # diameter
    @test_throws ArgumentError diameter(H, N(1 // 2))
    for res in ((@inferred diameter(H)), (@inferred diameter(H, Inf)))
        @test res === N(Inf)
    end
    @test_broken diameter(H, 2)  # TODO this should work

    # dim
    @test (@inferred dim(H)) == 2
    @test (@inferred dim(H3)) == 3

    # eltype
    @test (@inferred eltype(H)) == N
    @test (@inferred eltype(typeof(H))) == N

    # extrema
    @test extrema(H) == (N[-Inf, -Inf], N[Inf, Inf])
    @test_throws DimensionMismatch extrema(H, 3)
    @test extrema(H, 1) == (N(-Inf), N(Inf))

    # halfspace_left & halfspace_right
    H2 = halfspace_left(N[1, 1], N[2, 2])
    @test isidentical(H2, HalfSpace(N[1, -1], N(0)))
    H2 = halfspace_right(N[1, 1], N[2, 2])
    @test isidentical(H2, HalfSpace(N[-1, 1], N(0)))

    # high
    @test high(H) == N[Inf, Inf]
    @test_throws DimensionMismatch high(H, 3)
    @test high(H, 1) == N(Inf)

    # isbounded
    @test !(@inferred isbounded(H))

    # isboundedtype
    @test !(@inferred isboundedtype(typeof(H)))

    # iscomplement
    for (H1a, H1b, eq) in [(HalfSpace(N[-1], N(0)), HalfSpace(N[1 // 2], N(0)), true),
                         (H, HalfSpace(N[-2, -4], N(-6)), true),
                         (H, HalfSpace(N[-2, 4], N(-6)), false),
                         (H, HalfSpace(N[-2, -4], N(6)), false)]
        @test LazySets.iscomplement(H1a, H1b) == eq
    end

    # isconvex
    @test @inferred isconvex(H)

    # isconvextype
    @test @inferred isconvextype(typeof(H))

    # isempty
    @test !(@inferred isempty(H))
    @test_broken !(@inferred isempty(H, true))  # TODO make this type-stable (witness)
    res, w = isempty(H, true)
    @test !res && w isa Vector{N} && w ∈ H

    # isoperation
    @test !(@inferred isoperation(H))

    # isoperationtype
    @test !(@inferred isoperationtype(typeof(H)))

    # ispolyhedral
    @test @inferred ispolyhedral(H)

    # ispolyhedraltype
    @test @inferred ispolyhedraltype(typeof(H))

    # ispolytopic
    @test !(@inferred ispolytopic(H))

    # ispolytopictype
    @test !(@inferred ispolytopictype(typeof(H)))

    # isuniversal
    @test !(@inferred isuniversal(H))
    @test_broken @inferred isuniversal(H, true)  # TODO make this type-stable (witness)
    res, w = isuniversal(H, true)
    @test !res && w isa Vector{N} && w ∉ H

    # is_tighter_same_dir_2D
    H2a = HalfSpace(N[1, 0], N(1))
    H2b = HalfSpace(N[2, 0], N(2))
    H2c = HalfSpace(N[2, 0], N(1))
    @test LazySets.is_tighter_same_dir_2D(H2a, H2b) &&
          LazySets.is_tighter_same_dir_2D(H2c, H2b) &&
          LazySets.is_tighter_same_dir_2D(H2a, H2a)
    @test !LazySets.is_tighter_same_dir_2D(H2a, H2b; strict=true) &&
          LazySets.is_tighter_same_dir_2D(H2c, H2b; strict=true) &&
          !LazySets.is_tighter_same_dir_2D(H2a, H2a; strict=true)

    # low
    @test low(H) == N[-Inf, -Inf]
    @test_throws DimensionMismatch low(H, 3)
    @test low(H, 1) == N(-Inf)

    # norm
    @test_throws ArgumentError norm(H, N(1 // 2))
    for res in ((@inferred norm(H)), (@inferred norm(H, Inf)))
        @test res === N(Inf)
    end
    @test_broken norm(H, 2)  # TODO this should work

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        P = polyhedron(H)
        @test P isa Polyhedra.DefaultPolyhedron
        if N != Float32
            @test P isa Polyhedra.DefaultPolyhedron{N}
        end
    end

    # radius
    @test_throws ArgumentError radius(H, N(1 // 2))
    for res in ((@inferred radius(H)), (@inferred radius(H, Inf)))
        @test res === N(Inf)
    end
    @test_broken radius(H, 2)  # TODO this should work

    # rectify
    @test_broken rectify(H)  # TODO this should work

    # reflect
    @test_broken (@inferred reflect(H)) == HalfSpace(N[-1, -2], N(3))  # TODO specialize
    @test (@inferred reflect(H)) == HPolyhedron([HalfSpace(N[-1, -2], N(3))])

    # tosimplehrep
    @test (@inferred tosimplehrep(H)) == (N[1 2], N[3])

    # triangulate
    @static if isdefined(@__MODULE__, :MiniQhull)  # TODO this should work without MiniQhull (then remove the dependency)
        @test_throws ArgumentError triangulate(H)
    end

    # triangulate_faces
    @test_throws DimensionMismatch triangulate_faces(H)
    @test_throws ArgumentError triangulate_faces(H3)

    # vertices_list
    @test_throws ArgumentError vertices_list(H)

    # vertices
    @test_throws ArgumentError vertices(H)

    # volume
    @test_broken (@inferred volume(H)) == N(Inf)  # TODO this should work

    # affine_map
    @test_throws DimensionMismatch affine_map(ones(N, 2, 3), H, N[1, 1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 2), H, N[1])
    @static if isdefined(@__MODULE__, :CDDLib)
        H2 = affine_map(N[1 0; 0 2], H, N[1, 1])
        @test isequivalent(H2, HalfSpace(N[1, 1], N(5)))
        H2 = affine_map(ones(N, 2, 2), H, N[1, 1])
        # @test isequivalent(H2, Line(N[0, 0], N[1, 1]))  # TODO Rational
        H2 = affine_map(N[1 0; 0 2; 0 0], H, N[1, 1, 3])
        @test isequivalent(H2,
                           HPolyhedron([HalfSpace(N[1, 1, 0], N(5)), HalfSpace(N[0, 0, 1], N(3)),
                                        HalfSpace(N[0, 0, -1], N(-3))]))
    end

    # distance (between point and set)
    @test_throws DimensionMismatch distance(H, N[0])
    @test_throws ArgumentError distance(H, N[0, 0]; p=N(1 // 2))
    x = N[0, 0]
    xs = [N[0, 0], N[2, 3]]
    ys = [N[0, 0], N[1, 1]]  # closest points in H
    for (x, y) in zip(xs, ys)
        @test (@inferred distance(x, H)) == (@inferred distance(H, x)) == @inferred distance(x, y)
        for p in N[2]
            @test (@inferred distance(x, H; p=p)) == (@inferred distance(H, x; p=p)) ==
                  @inferred distance(x, y; p=p)
        end
        for p in N[1, Inf]
            @test_broken distance(x, H; p=p) == distance(H, x; p=p) == distance(x, y; p=p)
        end
    end

    # in (part 1)
    @test_throws DimensionMismatch N[0] ∈ H
    v1 = N[0, 0]
    @test @inferred v1 ∈ H
    v2 = N[2, 2]
    @test @inferred v2 ∉ H

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0], H)
    v3 = N[1, 1]
    @test_throws ArgumentError is_interior_point(v1, H; ε=N(0))
    @test_throws ArgumentError is_interior_point(v1, H; p=N(1 // 2))
    if N <: AbstractFloat
        @test_broken @inferred is_interior_point(v1, H)  # TODO specialize or make this type-stable (witness)
        @test is_interior_point(v1, H)
        @test_broken !(@inferred is_interior_point(v2, H))  # TODO specialize or make this type-stable (witness)
        @test !is_interior_point(v2, H)
        @test_broken !(@inferred is_interior_point(v3, H))  # TODO specialize or make this type-stable (witness)
        @test_broken !is_interior_point(v3, H)  # TODO this should work
    else
        @test_throws ArgumentError is_interior_point(v1, H)
        @test is_interior_point(v1, H; ε=1 // 100)
        @test !is_interior_point(v2, H; ε=1 // 100)
        @test !is_interior_point(v3, H; ε=1 // 100)
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([0.0, 0.0], H)
    end

    # linear_map (part 1)
    @test_throws DimensionMismatch linear_map(ones(N, 2, 1), H)
    H2 = linear_map(N[1 0; 0 2], H)
    @test isequivalent(H2, HalfSpace(N[1, 1], N(3)))
    @static if isdefined(@__MODULE__, :Polyhedra)
        X = linear_map(ones(N, 2, 2), H)
        @test isequivalent(X, Line(N[0, 0], N[1, 1]))
        X = linear_map(zeros(N, 2, 2), H)  # zero map
        @test isequivalent(X, Singleton(N[0, 0]))
    end
    # TODO add test for Universe result

    # linear_map_inverse
    @test_throws AssertionError LazySets.linear_map_inverse(ones(N, 1, 1), H)
    # invertible map
    H2 = LazySets.linear_map_inverse(N[1 0; 0 1], H)
    @test isequivalent(H2, H)
    # noninvertible map
    M = ones(N, 2, 1)
    X = LazySets.linear_map_inverse(M, H)
    A, b = tosimplehrep(H)
    Y = HPolyhedron(A * M, b)
    @test X isa LazySet{N} && isequivalent(X, Y)
    # TODO add tests for Universe and EmptySet result

    # permute
    @test_throws DimensionMismatch permute(H, [1])
    @test_throws DimensionMismatch permute(H, [1, -1])
    @test_throws DimensionMismatch permute(H, [1, 3])
    @test_throws ArgumentError permute(H, [1, 1])
    @test isidentical((@inferred permute(H, [1, 2])), H)
    H2 = @inferred permute(H, [2, 1])
    @test isidentical(H2, HalfSpace(N[2, 1], N(3)))

    # project
    @test_throws DimensionMismatch project(H, [1, 2, 3])
    @test_throws DimensionMismatch project(H, [1, -1])
    @test_throws DimensionMismatch project(H, [1, 3])
    @test_throws ArgumentError project(H, [1, 1])
    H2 = project(HalfSpace(N[0, 1], N(2)), [2])
    @test isidentical(H2, HalfSpace(N[1], N(2)))
    X = project(H, [2])
    @test X isa Universe{N} && dim(X) == 1

    # sample
    res = @inferred sample(H)
    @test res isa Vector{N} && res ∈ H
    res = @inferred sample(H, 2)
    @test res isa Vector{Vector{N}} && length(res) == 2 && all(x ∈ H for x in res)

    # scale
    H2 = scale(N(2), H)  # TODO specialize so that `@inferred` does not show Universe result
    @test isidentical(H2, HalfSpace(N[1//2, 1], N(3)))
    # degenerate case
    @static if isdefined(@__MODULE__, :Polyhedra)
        X = scale(N(0), H)
        @test isequivalent(X, ZeroSet{N}(2))
    end
    # scale!
    @test_broken @inferred scale!(N(2), H)  # TODO this should be supported
    # degenerate case
    # @test_throws ArgumentError scale!(N(0), H)

    # support_function
    @test_throws DimensionMismatch ρ(N[1], H)
    res = @inferred ρ(N[1, 1], H)
    @test res === N(Inf)
    res = @inferred ρ(N[1, 2], H)
    @test res === N(3)
    res = @inferred ρ(N[0, 0], H)
    @test res === N(0)

    # support_vector
    @test_throws DimensionMismatch σ(N[1], H)
    @test_throws ArgumentError σ(N[1, 1], H)
    res = @inferred σ(N[1, 2], H)
    @test res == N[3, 0]

    # translate
    @test_throws DimensionMismatch translate(H, N[1])
    H2 = @inferred translate(H, N[1, 2])
    @test isidentical(H2, HalfSpace(N[1, 2], N(8)))
    # translate!
    # @test_throws DimensionMismatch translate!(H, N[1])  # TODO this should be supported (but does it always work? we can only modify `H.a`)
    H2 = copy(H)
    @test_broken @inferred translate!(H2, N[1, 2])

    # cartesian_product
    H2 = @inferred cartesian_product(H, H1)
    @test isequivalent(H2, HPolyhedron([HalfSpace(N[0, 0, 3], N(4)), HalfSpace(N[1, 2, 0], N(3))]))
    H2 = @inferred cartesian_product(H1, H)
    @test isequivalent(H2, HPolyhedron([HalfSpace(N[0, 1, 2], N(3)), HalfSpace(N[3, 0, 0], N(4))]))

    # convex_hull (binary)
    @test_throws DimensionMismatch convex_hull(H, H3)
    @test_broken convex_hull(H, H)  # TODO this should work
    # H2 = convex_hull(H, H)
    # @test isidentical(H2, H)
    H2 = HalfSpace(N[0, 3], N(4))
    @test_broken convex_hull(H, H2)  # TODO this should work
    # for X in (convex_hull(H, H2), convex_hull(H2, H))
    #     X = convex_hull(H, H)
    #     @test X isa Universe{N} && dim(X) == 2
    # end

    # difference
    @test_throws DimensionMismatch difference(H, H3)
    # TODO better implementation: either H1 or empty set or H1 ∩ H2⁻¹
    # full set
    X = difference(H, complement(H))
    @test X isa LazySet{N} && isequivalent(X, H)
    # empty set
    X = difference(H, H)
    @test X == UnionSetArray{N,LazySet{N}}(LazySet{N}[])
    # one polyhedron
    H2a = HalfSpace(N[0, 3], N(4))
    X = difference(H, H2a)
    @test X isa LazySet{N} && isequivalent(X, HPolyhedron([H, complement(H2a)]))

    # distance (between two sets)
    @test_broken distance(H, H)  # TODO this should work

    # exact_sum
    @test_throws DimensionMismatch exact_sum(H, H3)
    @static if isdefined(@__MODULE__, :Polyhedra)
        H2 = HalfSpace(N[2, 3], N(4))
        # TODO finish this test (not tested yet)
        for H2b in ((@inferred exact_sum(H, H2)), @inferred exact_sum(H2, H))
            @test H2b isa Hyperrectangle{N} &&
                isequivalent(H2b, Hyperrectangle(N[5, -3], N[3, 6]))
        end
    end

    # intersection
    @test_throws DimensionMismatch intersection(H, H3)
    # disjoint
    X = intersection(H, HalfSpace(N[-1, -2], N(-4)))
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    # overlapping
    H2 = HalfSpace(N[-1, -2], N(-2))
    X = intersection(H, H2)
    @test isequivalent(X, HPolyhedron([H, H2]))
    # identical
    X = intersection(H, H)
    @test isequivalent(X, H)

    # isapprox
    @test @inferred H ≈ H
    res = (H ≈ translate(H, N[1 // 100000000, 0]))
    if N <: AbstractFloat
        @test res  # below default tolerance for AbstractFloat
    else
        @test !res  # zero default tolerance for Rational
    end
    @test !(@inferred H ≈ translate(H, N[1 // 100, 0]))  # above default tolerance for all types
    @test !(@inferred H ≈ H3) && !(@inferred H3 ≈ H) && !(@inferred H ≈ P) && !(@inferred P ≈ H)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(H, H3)
    # disjoint
    H2 = HalfSpace(N[-1, -2], N(-4))
    for (H2a, H2b) in ((H, H2), (H2, H))
        @test_broken @inferred isdisjoint(H2a, H2b)  # TODO make this type-stable (witness)
        @test isdisjoint(H2a, H2b)
        @test_broken @inferred isdisjoint(H2a, H2b, true)  # TODO make this type-stable (witness)
        res, w = isdisjoint(H2a, H2b, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    # overlapping
    H2 = HalfSpace(N[-1, -2], N(-2))
    for (H2a, H2b) in ((H, H), (H, H2), (H2, H))
        @test !isdisjoint(H2a, H2b)
        res, w = isdisjoint(H2a, H2b, true)
        @test !res && w isa Vector{N} && w ∈ H2a && w ∈ H2b
    end

    # isequal
    @test @inferred H == H
    @test (@inferred H != H3) && (@inferred H3 != H) && (@inferred H != P) && @inferred P != H

    # isequivalent
    @test_throws DimensionMismatch isequivalent(H, H3)
    @test_broken @inferred isequivalent(H, H)  # TODO make this type-stable (witness)
    @test isequivalent(H, H)
    @test !isequivalent(H, HalfSpace(N[-1, -2], N(-2)))
    @test isequivalent(H, P) && isequivalent(P, H)

    # isstrictsubset
    @test_throws DimensionMismatch H ⊂ H3
    # equivalent
    for X in (H, P)
        @test_broken !(@inferred X ⊂ H)  # TODO make this type-stable (witness)
        @test !(X ⊂ H)
        @test_broken @inferred ⊂(X, H, true)  # TODO make this type-stable (witness)
        res, w = ⊂(X, H, true)
        @test !res && w isa Vector{N} && isempty(w)
    end
    # no subset
    H2 = HalfSpace(N[-1, -2], N(-2))
    @test !(H2 ⊂ H)
    @test_broken ⊂(H2, H, true)  # TODO this should work
    # res, w = ⊂(H2, H, true)
    # @test !res && w isa Vector{N} && w ∈ H2 && w ∉ H
    # strict subset
    H2 = HalfSpace(N[1, 2], N(4))
    @test H ⊂ H2
    res, w = ⊂(H, H2, true)
    @test res && w isa Vector{N} && w ∉ H && w ∈ H2

    # issubset
    @test_throws DimensionMismatch H ⊆ H3
    for X in (H, P, H2 = HalfSpace(N[1, 2], N(4)))
        @test_broken @inferred H ⊆ X  # TODO make this type-stable (witness)
        @test H ⊆ X
        @test_broken @inferred ⊆(H, X, true)  # TODO make this type-stable (witness)
        res, w = ⊆(H, X, true)
        @test res && w isa Vector{N} && w == N[]
    end
    H2 = H2 = HalfSpace(N[1, 2], N(0))
    @test H ⊈ H2
    res, w = ⊆(H, H2, true)
    @test !res && w isa Vector{N} && w ∈ H && w ∉ H2

    # linear_combination
    @test_throws DimensionMismatch linear_combination(H, H3)
    @test_broken linear_combination(H, Xnc) isa LazySet{N}  # TODO implement `linear_combination` for non-convex sets
    @test_broken linear_combination(Xnc, H) isa LazySet{N}
    @test_broken linear_combination(H, H)  # TODO this should work
    # for X in (linear_combination(H, H), linear_combination(H, P), linear_combination(P, H))
    #     @test X isa LazySet{N} && isequivalent(X, H)
    # end
    # H2 = HalfSpace(N[1, -2], N(0))
    # for X in (linear_combination(H, H2), linear_combination(H2, H))
    #     @test X isa LazySet{N} && ...
    # end

    # minkowski_difference
    @test_throws DimensionMismatch minkowski_difference(H, H3)
    @test_broken minkowski_difference(H, H)  # TODO support this

    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(H, H3)
    H2 = HalfSpace(N[1, -2], N(0))
    @static if isdefined(@__MODULE__, :Polyhedra)
        # TODO finish this test (not tested yet)
        for H2b in ((@inferred minkowski_sum(H, H2)), @inferred minkowski_sum(H2, H))
            @test H2b isa Hyperrectangle{N} &&
                isequivalent(H2b, Hyperrectangle(N[5, -3], N[3, 6]))
        end
    end
end
    ### TODO old tests below

    # test concrete linear map of a half-space
    H = HalfSpace(N[1, -1], N(0)) # x <= y
    M = N[1 0; 0 0]  # non-invertible matrix
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test_throws ArgumentError linear_map(M, H, algorithm="vrep")
    end
    M = N[2 2; 0 1]  # invertible matrix
    @test linear_map(M, H) == HalfSpace(N[0.5, -2.0], N(0.0))
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :CDDLib)
        M = hcat(N[1 1])  # universal
        @test linear_map(M, H) == Universe{N}(1)
        M = zeros(N, 2, 2)  # result is a singleton
        X = linear_map(M, H)
        @test X isa HPolyhedron && isequivalent(X, ZeroSet{N}(2))

        if N isa AbstractFloat
            @test linear_map(N[1 1], H) == Universe{N}(1)
        end
    end

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
    # infeasible constraints
    clist = [HalfSpace(N[-1], N(-1)), HalfSpace(N[1], N(0))]
    clist2 = remove_redundant_constraints(clist)
    @test clist2 isa Vector{HalfSpace{N,Vector{N}}} && isempty(clist2)
end

for N in @tN([Float64, Float32])
    H = @inferred HalfSpace(N[1, 2], N(3))

    # normalize
    H2 = normalize(H)
    @test norm(H2.a) ≈ N(1) && H2.b == H.b / norm(H.a)
    @test normalize(H, N(1)) == HalfSpace(N[1 // 3, 2 // 3], N(1))
    @test normalize(H, N(Inf)) == HalfSpace(N[1 // 2, 1], N(3 // 2))

    # rand
    H2 = rand(HalfSpace; N=N)
    @test H2 isa HalfSpace{N} && dim(H2) == 2
    H2 = rand(HalfSpace; N=N, dim=3)
    @test H2 isa HalfSpace{N} && dim(H2) == 3

    # rationalize
    H2 = rationalize(H)
    @test isa(H2, HalfSpace{Rational{Int},Vector{Rational{Int}}})
    H2.a == Rational{Int}[1 // 1, 2 // 1]
    H2 = rationalize(BigInt, H)
    @test isa(H2, HalfSpace{Rational{BigInt},Vector{Rational{BigInt}}})
    H2.a == Rational{BigInt}[1 // 1, 2 // 1]

    # exponential_map
    @test_throws DimensionMismatch exponential_map(ones(N, 1, 1), H)
    @test_throws DimensionMismatch exponential_map(ones(N, 3, 2), H)
    @test_broken @inferred exponential_map(ones(N, 2, 2), H)  # TODO this should work because  the matrix is invertible
    H2 = exponential_map(N[1 0; 0 2], H)
    @test isidentical(H2, HalfSpace(N[1 / ℯ, 2 / ℯ^2], N(3)))

    # in (part 2)
    # check robustness (see #2312)
    x = N[0.07768723948819561, -0.5762273280928935, 0.28897399484750297, 1.9299362784322858]
    H2 = HalfSpace(N[-0.09291863543681655, -0.2176689899601838, -0.07453829739226348,
                     0.048948632014371496], N(0.1911363393469332))
    @test x ∈ H2

    # linear_map (part 2)
    X = linear_map(N[1 0; 0 2; 0 0], H)
    @test isequivalent(X,
                       HPolyhedron([HalfSpace(N[1, 1, 0], N(3)), HalfSpace(N[0, 0, 1], N(0)),
                                    HalfSpace(N[0, 0, -1], N(0))]))
end

for N in [Float64]
    # constructor
    @static if isdefined(@__MODULE__, :Symbolics)
        # 1 variable
        vars = @variables x
        @test HalfSpace(x <= 2.0, vars) == HalfSpace([1.0], 2.0)

        # 2 variables
        vars = @variables x y
        @test_throws ArgumentError HalfSpace(2x == 5y, vars)

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

        # because `get_vars` returns variables [y, x], both tests below require `vars` to pass
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

    # convert
    @static if isdefined(@__MODULE__, :SymEngine)
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

    # _ishalfspace
    @static if isdefined(@__MODULE__, :SymEngine)
        res = all(LazySetsSymEngineExt._ishalfspace.([:(x1 <= 0), :(x1 < 0), :(x1 > 0), :(x1 >= 0)]))
        res &= !LazySetsSymEngineExt._ishalfspace(:(x1 = 0))
        res &= LazySetsSymEngineExt._ishalfspace(:(2 * x1 <= 4))
        res &= LazySetsSymEngineExt._ishalfspace(:(6.1 <= 5.3 * f - 0.1 * g))
        res &= !LazySetsSymEngineExt._ishalfspace(:(2 * x1^2 <= 4))
        res &= !LazySetsSymEngineExt._ishalfspace(:(x1^2 > 4 * x2 - x3))
        res &= LazySetsSymEngineExt._ishalfspace(:(x1 > 4 * x2 - x3))
        @test res
    end
end
