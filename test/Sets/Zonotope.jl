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
    P = VPolygon([N[3, 4], N[5, 8], N[-1, 0], N[-3, -4]])  # equivalent set
    Xnc = UnionSet(P, BallInf(N[3, 0], N(1)))  # nonconvex set

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

    # complement
    if isdefined(@__MODULE__, :Polyhedra) || N <: AbstractFloat
        X = complement(Z3)
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
    Z2 = concretize(Z)
    @test isidentical(Z, Z2)

    # constrained_dimensions
    @test constrained_dimensions(Z) == 1:2

    if isdefined(@__MODULE__, :Polyhedra) || N <: AbstractFloat
        # constraints_list
        cs = constraints_list(Z3)
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
        res = constraints_list(Z2)
        @test res isa Vector{<:HalfSpace{N}} &&
              ispermutation(res, [HalfSpace(N[1], N(1)), HalfSpace(N[-1], N(1))])
        # flat zonotope (#3209)
        Z2 = Zonotope(zeros(N, 3), N[1 1 0; 0 1 1; 0 0 0])
        cs2 = constraints_list(Z2)
        clist = [HalfSpace(N[1, 0, 0], N(2)), HalfSpace(N[-1, 0, 0], N(2)),
                 normalize(HalfSpace(N[1, -1, 0], N(2))), normalize(HalfSpace(N[-1, 1, 0], N(2))),
                 HalfSpace(N[0, -1, 0], N(2)), HalfSpace(N[0, 1, 0], N(2)),
                 HalfSpace(N[0, 0, 1], N(0)), HalfSpace(N[0, 0, -1], N(0))]
        @test clist isa Vector{<:HalfSpace{N}}
        if N <: AbstractFloat
            @test ispermutation(cs2, clist)
        else
            @test all(isequivalent.(cs2, clist))  # rely on order of constraints for this test
        end

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

    # triangulate_faces
    @test_throws DimensionMismatch triangulate_faces(Z)
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :GeometryBasics)
        res = triangulate_faces(Z3)
        @test res isa Tuple{Matrix{Float32},Vector{Tuple{Int,Int,Int}}}
        @test length(res[2]) == 12
        # no precise test here; the `triangulate_faces` implementation is tested elsewhere
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

    # affine_map
    @test_throws DimensionMismatch affine_map(ones(N, 2, 1), Z, N[1, 1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 2), Z, N[1])
    Z2 = affine_map(ones(N, 1, 2), Z, N[1])
    @test isidentical(Z2, Zonotope(N[4], hcat(N(10))))
    Z2 = affine_map(ones(N, 3, 2), Z, N[1, 2, 3])
    @test isidentical(Z2, Zonotope(N[4, 5, 6], N[3 7; 3 7; 3 7]))

    # distance (between point and set)
    @test_broken distance(Z, N[0, 0])  # TODO this should work
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

    # exponential_map
    @test_throws DimensionMismatch exponential_map(ones(N, 1, 1), Z)
    @test_throws DimensionMismatch exponential_map(ones(N, 2, 1), Z)

    # in
    @test_throws DimensionMismatch N[0] ∈ Z
    @test N[1, 2] ∈ Z && N[5, 8] ∈ Z && N[3, 3] ∉ Z

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0], Z)
    @test_throws ArgumentError is_interior_point(N[0, 0], Z; ε=N(0))
    @test_throws ArgumentError is_interior_point(N[0, 0], Z; p=N(1 // 2))
    if N <: AbstractFloat
        @test is_interior_point(N[1, 2], Z)
        @test !is_interior_point(N[5, 8], Z)
        @test !is_interior_point(N[3, 3], Z)
    else
        @test_throws ArgumentError is_interior_point(N[0, 0], Z)
        @test is_interior_point(N[1, 2], Z; ε=1 // 100)
        @test !is_interior_point(N[5, 8], Z; ε=1 // 100)
        @test !is_interior_point(N[3, 3], Z; ε=1 // 100)
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([1.0], Z)
    end

    # linear_map
    @test_throws DimensionMismatch linear_map(ones(N, 1, 1), Z)
    Z2 = linear_map(N[1 0; 0 2], Z)
    @test isidentical(Z2, Zonotope(N[1, 4], N[1 3; 4 8]))
    Z2 = linear_map(zeros(N, 1, 2), Z)  # zero map
    @test isidentical(Z2, Zonotope(N[0], zeros(N, 1, 0)))
    # 1D simplifies to 1 generator (except for the zero map)
    Z2 = linear_map(N[-1 1;], Z)
    @test isidentical(Z2, Zonotope(N[1], hcat(N[2])))
    # linear_map!
    @test_broken linear_map!(N[1 0; 0 2], Z2) isa Zonotope  # TODO this should work
    # Z2 = copy(Z)
    # linear_map!(N[1 0; 0 2], Z2)
    # @test isidentical(Z2, Zonotope(N[1, 4], N[1 3; 4 8]))

    # linear_map_inverse
    @test_throws AssertionError LazySets.linear_map_inverse(ones(N, 1, 1), Z)
    X = LazySets.linear_map_inverse(N[1 0; 0 1], Z)  # TODO this could be a Zonotope
    @test_broken isidentical(X, Z)
    @test X isa LazySet{N} && isequivalent(X, Z)
    M = ones(N, 2, 1)
    X = LazySets.linear_map_inverse(M, Z)  # TODO this could be a Zonotope
    @test_broken isidentical(X, Zonotope(N[-1], hcat(N[1])))
    A, b = tosimplehrep(Z)
    Y = HPolytope(A * M, b)
    @test X isa LazySet{N} && isequivalent(X, Y)

    # permute
    @test_throws DimensionMismatch permute(Z, [1])
    @test_throws DimensionMismatch permute(Z, [1, -1])
    @test_throws DimensionMismatch permute(Z, [1, 3])
    @test_throws ArgumentError permute(Z, [1, 1])
    Z2 = permute(Z, [1, 2])
    @test isidentical(Z2, Z)
    Z2 = permute(Z, [2, 1])
    @test isidentical(Z2, Zonotope(N[2, 1], N[2 4; 1 3]))

    # project
    @test_throws DimensionMismatch project(Z, [1, 2, 3])
    @test_throws DimensionMismatch project(Z, [1, -1])
    @test_throws DimensionMismatch project(Z, [1, 3])
    @test_throws ArgumentError project(Z, [1, 1])
    Z2 = project(Z, [1])
    @test isidentical(Z2, Zonotope(N[1], hcat(N[4])))
    # 1D with zero generators (#2147)
    Z2 = project(convert(Zonotope, BallInf(N[0, 0], N(1))), [1])
    @test isidentical(Z2, Zonotope(N[0], hcat(N[1])))
    # removing zero generators in projection can be deactivated
    Z2a = Zonotope(zeros(N, 3), N[1 0 3; 4 0 6; 0 0 0])
    Z2 = project(Z2a, [1, 2])
    @test isidentical(Z2, Zonotope(zeros(N, 2), N[1 3; 4 6]))
    Z2 = project(Z2a, [1, 2]; remove_zero_generators=false)
    @test isidentical(Z2, Zonotope(zeros(N, 2), N[1 0 3; 4 0 6]))

    # sample
    if N <: AbstractFloat
        res = sample(Z)
        @test res isa Vector{N} && res ∈ Z
        res = sample(Z, 2)
        @test res isa Vector{Vector{N}} && length(res) == 2 && all(x ∈ Z for x in res)
    else
        @test_broken sample(Z) isa Vector{N}  # TODO this should work
        @test_broken sample(Z, 2) isa Vector{Vector{N}}
    end

    # scale
    Z2 = scale(N(2), Z)
    @test isidentical(Z2, Zonotope(N[2, 4], N[2 6; 4 8]))
    # degenerate case
    Z2 = scale(N(0), Z)
    @test isidentical(Z2, Zonotope(N[0, 0], zeros(N, 2, 2)))  # TODO remove zero generators
    # scale!
    Z2 = copy(Z)
    scale!(N(2), Z2)
    @test isidentical(Z2, Zonotope(N[2, 4], N[2 6; 4 8]))

    # split
    # split into two zonotopes
    @test_throws AssertionError split(Z, 0)
    Z2a, Z2b = split(Z, 2)
    @test isidentical(Z2a, Zonotope(N[-1 // 2, 0], N[1 3//2; 2 2]))
    @test isidentical(Z2b, Zonotope(N[5 // 2, 4], N[1 3//2; 2 2]))
    # split into multiple zonotopes
    Z2a = Zonotope(N[0, 0], N[1 1; -1 1])
    Z2a, Z2b, Z2c, Z2d = split(Z2a, [1, 2], [1, 1])
    G = N[1//2 1//2; -1//2 1//2]
    ispermutation([Z2a, Z2b, Z2c, Z2d],
                  [Zonotope(N[-1, 0], G), Zonotope(N[0, 1], G),
                   Zonotope(N[0, -1], G), Zonotope(N[1, 0], G)])
    @static if isdefined(@__MODULE__, :StaticArrays)
        using StaticArrays: SVector, SMatrix

        Z2 = Zonotope(SVector{2}(N[0, 0]), SMatrix{2,2}(N[1 1; -1 1]))
        Z2a, Z2b = split(Z2, 1)
        @test Z2a ⊆ Z2 && Z2b ⊆ Z2
    end

    # support_function
    @test_throws DimensionMismatch ρ(N[1], Z)
    res = ρ(N[1, 1], Z)
    @test res isa N && res == N(13)
    res = ρ(N[-1, -1], Z)
    @test res isa N && res == N(7)

    # support_vector
    @test_throws DimensionMismatch σ(N[1], Z)
    res = σ(N[1, 1], Z)
    @test res isa Vector{N} && res == N[5, 8]
    res = σ(N[-1, -1], Z)
    @test res isa Vector{N} && res == N[-3, -4]

    # translate
    @test_throws DimensionMismatch translate(Z, N[1])
    Z2 = translate(Z, N[1, 2])
    @test isidentical(Z2, Zonotope(N[2, 4], N[1 3; 2 4]))
    # translate!
    @test_throws DimensionMismatch translate!(Z, N[1])
    Z2 = copy(Z)
    translate!(Z2, N[1, 2])
    @test isidentical(Z2, Zonotope(N[2, 4], N[1 3; 2 4]))

    # cartesian_product
    Z2 = Zonotope(N[3], hcat(N[1]))
    Z2b = cartesian_product(Z, Z2)
    @test isidentical(Z2b, Zonotope(N[1, 2, 3], N[1 3 0; 2 4 0; 0 0 1]))
    Z2b = cartesian_product(Z2, Z)
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
    @test_broken difference(Z, Z3) isa LazySet{N}  # TODO this should work (add more tests later)

    # distance (between two sets)
    @test_broken distance(Z, Z3) isa LazySet{N}  # TODO this should work (add more tests later)

    # exact_sum
    @test_throws DimensionMismatch exact_sum(Z, Z3)
    Z2 = Zonotope(N[4, -2], N[1 0; 0 2])
    for Z2b in (exact_sum(Z, Z2), exact_sum(Z2, Z))
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
    @test X isa LazySet{N} && isequivalent(X, Y)
    # intersection with HalfSpace
    H = HalfSpace(N[1, 0], N(-4))  # disjoint
    X = intersection(H, Z)
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    H = HalfSpace(N[1, 0], N(6))  # fully contained
    X = intersection(H, Z)
    @test isidentical(X, Z)
    H = HalfSpace(N[1, 0], N(0))  # overlapping
    X = VPolygon([N[-3, -4], N[0, 0], N[0, 4 // 3], N[-1, 0]])
    @test isequivalent(intersection(H, Z), X)

    # isapprox
    @test Z ≈ Z
    res = (Z ≈ translate(Z, N[1 // 100000000, 0]))
    if N <: AbstractFloat
        @test res  # below default tolerance for AbstractFloat
    else
        @test !res  # zero default tolerance for Rational
    end
    @test !(Z ≈ translate(Z, N[1 // 1000, 0]))  # above default tolerance for all types
    @test !(Z ≈ Z3) && !(Z3 ≈ Z) && !(Z ≈ P) && !(P ≈ Z)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(Z, Z3)
    # disjoint
    Z2 = Zonotope(N[2, -2], N[1 0; 0 1])
    @test isdisjoint(Z, Z2) && isdisjoint(Z2, Z)
    for (pair, Z) in ((isdisjoint(Z, Z2, true), Z2), (isdisjoint(Z2, Z, true), Z2))
        res, w = pair
        @test res && w isa Vector{N} && isempty(w)
    end
    # overlapping
    Z2 = Zonotope(N[2, 2], N[1 0; 0 1])
    @test !isdisjoint(Z, Z) && !isdisjoint(Z, Z2) && !isdisjoint(Z2, Z)
    @test_throws ErrorException isdisjoint(Z, Z, true)  # TODO this should be an ArgumentError
    # for (pair, Z) in ((isdisjoint(Z, Z, true), Z), (isdisjoint(Z, Z2, true), Z2),  # TODO this should work
    #                   (isdisjoint(Z2, Z, true), Z2))
    #     res, w = pair
    #     @test !res && w isa Vector{N} && w ∈ Z && w ∈ Z
    # end
    # tolerance
    if N == Float64
        Z2 = Zonotope(N[6 + 1e-9, 8], N[1 0; 0 1])
        @test !isdisjoint(Z, Z2)
        LazySets.set_rtol(Float64, 1e-10)
        @test_broken isdisjoint(Z, Z2)  # TODO this should work
        # restore tolerance
        LazySets.set_rtol(Float64, LazySets.default_tolerance(Float64).rtol)
    end
    # TODO revise
    # OLD TESTS
    Z1 = Zonotope(N[1, 1], N[1 1; -1 1])
    Z2 = Zonotope(N[-2, -1], Matrix{N}(I, 2, 2))
    # isdisjoint with a hyperplane
    H1 = Hyperplane(N[1, 1], N(3))
    intersection_empty, point = isdisjoint(Z1, H1, true)
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
    # intersection with a hyperplane
    Z1 = Zonotope(N[1, 1], N[1 1; -1 1])
    H1 = Hyperplane(N[1, 1], N(3))
    intersection_empty, point = isdisjoint(Z1, H1, true)
    @test !isdisjoint(Z1, H1) && !intersection_empty
    H2 = Hyperplane(N[1, 1], N(-11))
    @test isdisjoint(Z1, H2) && isdisjoint(Z1, H2, true)[1]
    @test !isdisjoint(H1, Z1)
    @test isdisjoint(H2, Z1) && isdisjoint(H2, Z1, true)[1]

    # isequal
    @test Z == Z
    @test Z != Z3 && Z3 != Z && Z != P && P != Z

    # isequivalent
    @test_throws DimensionMismatch isequivalent(Z, Z3)
    @test isequivalent(Z, Z)
    @test !isequivalent(Z, Zonotope(N[1, 2], N[1 0; 0 1]))
    @test isequivalent(Z, P) && isequivalent(P, Z)

    # isstrictsubset
    @test_throws DimensionMismatch Z ⊂ Z3
    @test_throws DimensionMismatch Z3 ⊂ Z
    for Z2 in (Z, P)
        @test !(Z2 ⊂ Z)
        res, w = ⊂(Z2, Z, true)
        @test !res && w isa Vector{N} && isempty(w)
    end
    Z2 = Zonotope(N[1, 2], N[1 0; 0 1])
    @test_broken ⊂(Z2, Z, true) isa Tuple  # TODO support witness production, then integrate in loop above
    Z2 = Zonotope(N[1, 2], N[4 0; 0 6])
    @test Z ⊂ Z2
    @test_broken ⊂(Z, Z2, true) isa Tuple  # TODO support witness production
    # res, w = ⊂(Z, Z2, true)
    # @test res && w isa Vector{N} && w ∉ Z && w ∈ Z2

    # TODO revise
    # issubset
    @test_throws DimensionMismatch Z ⊆ Z2
    @test_throws DimensionMismatch Z2 ⊆ Z
    for Z2 in (Z, B)
        @test Z ⊆ Z2
        res, w = ⊆(Z, Z2, true)
        @test res && w isa Vector{N} && w == N[]
    end
    for Z2 in (Interval(N(0), N(1)), Interval(N(1), N(3)))
        @test Z ⊈ Z2
        res, w = ⊆(Z, Z2, true)
        @test !res && w isa Vector{N} && w ∈ Z && w ∉ Z2
    end
    # OLD TESTS
    Z = Zonotope(N[0, 0], N[1 1; -1 1])
    H1 = Hyperrectangle(; low=N[-2, -2], high=N[2, 2])
    H2 = Hyperrectangle(; low=N[-2, -2], high=N[2, 0])
    @test issubset(Z, H1)
    @test !issubset(Z, H2)

    # linear_combination
    @test_throws DimensionMismatch linear_combination(Z, Z3)
    @test_broken linear_combination(Z, Xnc)
    @test_broken linear_combination(Xnc, Z)
    for X in (linear_combination(Z, Z), linear_combination(Z, P), linear_combination(P, Z))
        @test X isa LazySet{N} && isequivalent(X, Z)
    end
    Z2 = Zonotope(N[4, -2], N[1 0; 0 2])
    Y = VPolygon([N[-3, -4], N[5, -4], N[5, 8], N[-1, 0]])
    for X in (linear_combination(Z, Z2), linear_combination(Z2, Z))
        @test X isa LazySet{N} && isequivalent(X, Y)
    end

    # minkowski_difference
    @test_throws DimensionMismatch minkowski_difference(Z, Z3)
    @test_throws DimensionMismatch minkowski_difference(Z3, Z)
    # empty difference
    Z2 = Zonotope(N[1, 2], N[4 0; 0 6])
    X = minkowski_difference(Z, Z2)
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    # nonempty difference
    Z2 = minkowski_difference(Z, Z)
    @test isidentical(Z2, Zonotope(N[0, 0], zeros(N, 2, 0)))
    # TODO revise below
    Z2 = Interval(N(1), N(3))
    Z = minkowski_difference(Z, Z2)
    @test isidentical(Z, Interval(N(-1), N(-1)))
    for Z2 in (minkowski_difference(Z, B), minkowski_difference(B, Z))
        @test isequivalent(Z2, Interval(N(0), N(0)))
    end
    # Universe
    U = Universe{N}(1)
    U2 = minkowski_difference(U, Z)
    @test U2 isa Universe{N} && dim(U2) == 1

    # TODO revise
    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(Z, Z2)
    @test_throws DimensionMismatch minkowski_sum(Z2, Z)
    # Interval + Interval = Interval
    Z2 = minkowski_sum(Z, Z)
    Z = Interval(N(0), N(4))
    @test isidentical(Z2, Z)
    # general
    for Z2 in (minkowski_sum(Z, B), minkowski_sum(B, Z))
        @test Z2 isa LazySet{N} && isequivalent(Z2, Z)
    end
    # OLD TESTS
    Z2 = Zonotope(N[-1, 1], Matrix{N}(I, 2, 2))
    Z3 = minkowski_sum(Z1, Z2)
    @test Z3.center == N[0, 2]
    @test Z3.generators == N[1 1 1 0; -1 1 0 1]
end

for N in @tN([Float64, Float32, Rational{Int}])
    # TODO revise this loop and merge with above loop

    # zero column in generators
    g = zeros(N, 2, 5)
    g[:, 3] = ones(N, 2)
    g[1, 2] = N(2)
    z = Zonotope(N[1, 2], g)
    @test size(z.generators) == (2, 5)
    zred = remove_zero_generators(z)
    @test size(zred.generators) == (2, 2)

    # reduce_order
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

        Z = Zonotope(N[2, 1], N[-1//2 3//2 1//2 1 0 1; 1//2 3//2 1 1//2 1 0])
        for method in [LazySets.COMB03(), LazySets.GIR05()]
            Zs = Zonotope(SVector{2}(Z.center), SMatrix{2,6}(Z.generators))
            @test reduce_order(Zs, 2, method) isa Zonotope{N,SVector{2,N},SMatrix{2,4,N,8}}
        end
    end

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
        @test isequivalent(Z, H)
    end

    # same for CartesianProductArray
    Z2 = convert(Zonotope, CartesianProductArray([h1, h2]))
    @test Z == Z2

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
end

for N in @tN([Float64, Float32])
    Z = Zonotope(N[1, 2], N[1 3; 2 4])

    # rand
    Z2 = rand(Zonotope; N=N)
    @test Z2 isa Zonotope{N} && dim(Z2) == 2
    Z2 = rand(Zonotope; N=N, dim=3)
    @test Z2 isa Zonotope{N} && dim(Z2) == 3
    Z2 = rand(Zonotope; N=N, num_generators=5)
    @test Z2 isa Zonotope{N} && dim(Z2) == 2 && ngens(Z2) == 5

    # rationalize
    Z2 = rationalize(Z)
    T = Rational{Int64}
    @test Z2 isa Zonotope{T} && isidentical(Z2, Zonotope(T[1, 2], T[1 3; 2 4]))

    # exponential_map
    Z2 = exponential_map(N[1 0; 0 2], Z)
    @test isidentical(Z2, Zonotope(N[exp(1), 2 * exp(2)], N[exp(1) 3*exp(1); 2*exp(2) 4*exp(2)]))
end

for N in @tN([Float64, Rational{Int}])
    Z = Zonotope(N[1, 2], N[1 3; 2 4])

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        X = polyhedron(Z)
        @test X isa Polyhedra.DefaultPolyhedron{N}
        @test ispermutation(constraints_list(convert(HPolytope, X)), constraints_list(Z))
    end

    # triangulate
    @static if isdefined(@__MODULE__, :MiniQhull)
        X = triangulate(Z)
        @test X isa UnionSetArray{N,<:VPolytope{N}} && length(X) == 2
        for Xi in X
            @test ispermutation(vertices_list(Xi), [N[3, 4], N[-1, 0], N[-3, -4]]) ||
                  ispermutation(vertices_list(Xi), [N[3, 4], N[-1, 0], N[5, 8]])
        end
    end

    # volume
    @static if isdefined(@__MODULE__, :Polyhedra)
        res = volume(Z)
        @test res isa N && res ≈ N(8)  # TODO this should work without Polyhedra in 2D
    end
end

for N in [Float64]
    Z3 = Zonotope(N[1, 1, 1], N[1 0 0; 0 2 0; 0 0 3])  # 3D zonotope, a hyperrectangle

    # chebyshev_center_radius
    @static if isdefined(@__MODULE__, :Polyhedra)
        c, r = chebyshev_center_radius(Z3)
        @test c isa Vector{N} && c == Z3.center
        @test r isa N && r == N(1)
    end
end

for N in [Float64]
    # TODO revise this loop and merge with above loop

    # conversion to HPolytope
    # 1D
    Z = Zonotope(N[0], Matrix{N}(I, 1, 1))
    X = HPolytope(constraints_list(Z))
    for d in [N[1], N[-1]]
        @test ρ(d, X) == ρ(d, Z)
    end
    # 2D
    Z = Zonotope(N[0, 0], Matrix{N}(I, 2, 2))
    X = HPolytope(constraints_list(Z))
    for d in BoxDiagDirections{N}(2)
        @test ρ(d, X) == ρ(d, Z)
    end
end
