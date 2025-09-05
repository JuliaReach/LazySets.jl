using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

function isidentical(::Universe, ::Universe)
    return false
end

function isidentical(U1::Universe{N}, U2::Universe{N}) where {N}
    return U1.dim == U2.dim
end

let
    # default Float64 constructor
    U = Universe(2)
    @test U isa Universe{Float64}
    @test U.dim == 2
end

for N in @tN([Float64, Float32, Rational{Int}])
    # auxiliary sets
    B = BallInf(ones(N, 2), N(1))
    E = EmptySet{N}(2)
    Z = ZeroSet{N}(2)
    Pnc = Polygon([N[0, 0], N[3, 0], N[1, 1], N[0, 3]])  # nonconvex
    Pe = HPolygon([HalfSpace(N[1, 0], N(0)), HalfSpace(N[-1, 0], N(-1)),  # empty
                   HalfSpace(N[0, 1], N(0)), HalfSpace(N[0, -1], N(0))])

    # constructor
    U = Universe{N}(2)
    @test U isa Universe{N}
    @test U.dim == 2
    U3 = Universe{N}(3)
    @test U3 isa Universe{N}
    @test U3.dim == 3

    # an_element
    x = an_element(U)
    @test x isa Vector{N} && length(x) == 2

    # area
    @test_throws DimensionMismatch area(Universe{N}(1))
    for X in (U, U3)
        res = area(X)
        @test res isa N && res == N(Inf)
    end

    # chebyshev_center_radius
    @static if isdefined(@__MODULE__, :Polyhedra)
        # behavior differs for Rational solver:
        # - prints `glp_exact: problem has no rows/columns`
        # - different error message (reported in Polyhedra.jl#352)
        @test_throws ErrorException chebyshev_center_radius(U)
    end

    # complement
    E2 = complement(U)
    @test E2 isa EmptySet{N} && dim(E2) == 2

    # concretize
    U2 = concretize(U)
    @test isidentical(U2, U)

    # constrained_dimensions
    x = constrained_dimensions(U)
    @test x isa Vector{Int} && isempty(x)

    # constraints_list
    clist = constraints_list(U)
    @test isempty(clist) && eltype(clist) <: HalfSpace{N}

    # constraints
    citer = constraints(U)
    @test isempty(citer) && eltype(citer) <: AbstractVector{N}

    # convex_hull (unary)
    U2 = convex_hull(U)
    @test isidentical(U2, U)

    # copy
    U2 = copy(U)
    @test isidentical(U2, U)

    # diameter
    @test_throws ArgumentError diameter(U, N(1 // 2))
    @test_throws ArgumentError diameter(U)
    @test_throws ArgumentError diameter(U, Inf)
    @test_throws ArgumentError diameter(U, 2)

    # dim
    @test dim(U) == 2
    @test dim(U3) == 3

    # eltype
    @test eltype(U) == N
    @test eltype(typeof(U)) == N

    # extrema
    @test extrema(U) == (N[-Inf, -Inf], N[Inf, Inf])
    @test_throws DimensionMismatch extrema(U, 3)
    @test extrema(U, 1) == (N(-Inf), N(Inf))

    # high
    @test high(U) == N[Inf, Inf]
    @test_throws DimensionMismatch high(U, 3)
    @test high(U, 1) == N(Inf)

    # isbounded
    @test !isbounded(U)

    # isboundedtype
    @test !isboundedtype(typeof(U))

    # isconvextype
    @test isconvextype(typeof(U))

    # isempty
    @test !isempty(U)
    res, w = isempty(U, true)
    @test !res && w isa Vector{N} && w ∈ U

    # isoperation
    @test !isoperation(U)

    # isoperationtype
    @test !isoperationtype(typeof(U))

    # ispolyhedral
    @test ispolyhedral(U)

    # ispolyhedraltype
    @test ispolyhedraltype(typeof(U))

    # ispolytopic
    @test !ispolytopic(U)

    # isuniversal
    @test isuniversal(U)
    res, w = isuniversal(U, true)
    @test res && w isa Vector{N} && isempty(w)

    # low
    @test low(U) == N[-Inf, -Inf]
    @test_throws DimensionMismatch low(U, 3)
    @test low(U, 1) == N(-Inf)

    # norm
    @test_throws ArgumentError norm(U, N(1 // 2))
    @test_throws ArgumentError norm(U)
    @test_throws ArgumentError norm(U, Inf)
    @test_throws ArgumentError norm(U, 2)

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        P = polyhedron(U)
        @test P isa Polyhedra.DefaultPolyhedron
        if N != Float32
            @test P isa Polyhedra.DefaultPolyhedron{N}
        end
        @test size(P.hrep.A) == (0, 2)
    end

    # radius
    @test_throws ArgumentError radius(U, N(1 // 2))
    @test_throws ArgumentError radius(U)
    @test_throws ArgumentError radius(U, Inf)
    @test_throws ArgumentError radius(U, 2)

    # rand
    @test rand(Universe; N=N) isa Universe{N}
    U2 = rand(Universe; N=N, dim=3)
    @test isidentical(U2, U3)

    # rectify
    P = rectify(U)
    Q = HPolyhedron([HalfSpace(N[-1, 0], N(0)), HalfSpace(N[0, -1], N(0))])
    @test P isa HPolyhedron{N} && isequivalent(P, Q)

    # reflect
    @test reflect(U) == U

    # singleton_list
    @test_throws ArgumentError singleton_list(U)

    # tosimplehrep
    C, d = tosimplehrep(U)
    @test C isa Matrix{N} && d isa Vector{N} && size(C) == (0, 2) && isempty(d)

    # triangulate
    @test_throws ArgumentError triangulate(U)

    # triangulate_faces
    @test_throws DimensionMismatch triangulate_faces(U)
    @test_throws ArgumentError triangulate_faces(U3)

    # vertices_list
    @test_throws ArgumentError vertices_list(U)

    # vertices
    @test_throws ArgumentError vertices(U)

    # volume
    x = volume(U)
    @test x isa N && x == N(Inf)

    # affine_map
    @test_throws DimensionMismatch affine_map(ones(N, 2, 3), U, N[1, 1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 2), U, N[1])
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :CDDLib)
        # TODO this should work, even without Polyhedra
        @test_broken affine_map(ones(N, 2, 2), U, N[1, 1])
        # U2 = affine_map(ones(N, 2, 2), U, N[1, 1])
        # @test isidentical(U2, U)
        # U2 = affine_map(ones(N, 3, 2), U, N[1, 1, 3])
        # @test isidentical(U2, U3)
    end

    # distance (between point and set)
    @test_throws DimensionMismatch distance(U, N[0])
    @test_throws ArgumentError distance(U, N[0, 0]; p=N(1 // 2))
    x = N[0, 0]
    for res in (distance(U, x), distance(x, U))
        @test res isa N && res == N(0)
    end

    # exponential_map
    @test_throws DimensionMismatch exponential_map(ones(N, 1, 1), U)
    @test_throws DimensionMismatch exponential_map(ones(N, 3, 2), U)

    # in
    @test_throws DimensionMismatch N[0] ∈ U
    @test N[0, 0] ∈ U

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0], U)
    @test_throws ArgumentError is_interior_point(N[0, 0], U; ε=N(0))
    @test_throws ArgumentError is_interior_point(N[0, 0], U; p=N(1 // 2))
    if N <: AbstractFloat
        @test is_interior_point(N[0, 0], U)
    else
        @test_throws ArgumentError is_interior_point(N[0, 0], U)
        @test is_interior_point(N[0, 0], U; ε=1 // 100)
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([0.0, 0.0], U)
    end

    # linear_map
    @test_throws DimensionMismatch linear_map(ones(N, 2, 1), U)
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :CDDLib)
        # TODO this should work, even without Polyhedra
        U2 = linear_map(N[1 0; 0 1], U)
        @test U2 isa HPolyhedron{N}  # TODO this should change
        @test_broken isidentical(U2, U)
        if VERSION < v"1.12"
            # TODO these should work with older versions, see below
            @test_broken linear_map(ones(N, 2, 2), U)
            @test_broken linear_map(zeros(N, 2, 2), U)
            @test_broken linear_map(N[1 0; 0 1; 0 0], U)
        else
            U2 = linear_map(ones(N, 2, 2), U)
            @test U2 isa HPolyhedron{N}
            @test isequivalent(U2, Line(N[0, 0], N[1, 1]))
            X = linear_map(zeros(N, 2, 2), U)  # zero map
            @test isequivalent(X, ZeroSet{N}(2))
            U2 = linear_map(N[1 0; 0 1; 0 0], U)
            @test U2 isa HPolyhedron{N}
            @test isequivalent(U2, Hyperplane(N[0, 0, 1], N(0)))
        end
    end

    # linear_map_inverse
    U2 = LazySets.linear_map_inverse(ones(N, 2, 3), U)
    @test isidentical(U2, U3)

    # permute
    @test_throws DimensionMismatch permute(U, [1])
    @test_throws DimensionMismatch permute(U, [1, -1])
    @test_throws DimensionMismatch permute(U, [1, 3])
    @test_throws ArgumentError permute(U, [1, 1])
    U2 = permute(U, [1, 2])
    @test isidentical(U2, U)
    U2 = permute(U, [2, 1])
    @test isidentical(U2, U)

    # project
    @test_throws DimensionMismatch project(U, [1, 2, 3])
    @test_throws DimensionMismatch project(U, [1, -1])
    @test_throws DimensionMismatch project(U, [1, 3])
    @test_throws ArgumentError project(U, [1, 1])
    U2 = project(U, [2])
    @test U2 isa Universe{N} && dim(U2) == 1

    # sample
    x = sample(U)
    @test x isa Vector{N} && length(x) == 2
    xs = sample(U, 2)
    @test xs isa Vector{Vector{N}} && length(xs) == 2
    for x in xs
        @test x isa Vector{N} && length(x) == 2
    end

    # scale
    U2 = scale(N(2), U)
    @test isidentical(U2, U)
    Z = scale(N(0), U)
    @test Z isa ZeroSet{N} && Z == ZeroSet{N}(2)
    # scale!
    U2 = copy(U)
    scale!(N(2), U2)
    @test isidentical(U2, U)
    @test_throws ArgumentError scale!(N(0), U2)

    # support_function
    @test_throws DimensionMismatch ρ(N[1], U)
    v = ρ(N[-1, 2], U)
    @test v isa N && v == N(Inf)
    v = ρ(N[2, 0], U)
    @test v isa N && v == N(Inf)
    @test ρ(N[0, 0], U) == N(0)

    # support_vector
    @test_throws DimensionMismatch σ(N[1], U)
    x = σ(N[-1, 2], U)
    @test x isa Vector{N} && x == N[-Inf, Inf]
    x = σ(N[2, 0], U)
    @test x isa Vector{N} && x == N[Inf, 0]
    x = σ(N[0, 0], U)
    @test x isa Vector{N} && x == N[0, 0]

    # translate
    @test_throws DimensionMismatch translate(U, N[1])
    U2 = translate(U, N[1, 2])
    @test isidentical(U2, U)
    # translate!
    @test_throws DimensionMismatch translate!(U, N[1])
    U2 = copy(U)
    translate!(U2, N[1, 2])
    @test isidentical(U2, U)

    # cartesian_product
    for U2 in (cartesian_product(U, U3), cartesian_product(U3, U))
        @test U2 isa Universe{N} && dim(U2) == 5
    end

    # convex_hull (binary)
    @test_throws DimensionMismatch convex_hull(U, U3)
    U2 = convex_hull(U, U)
    @test isidentical(U2, U)
    for U2 in (convex_hull(U, Pnc), convex_hull(Pnc, U), convex_hull(U, E), convex_hull(E, U))
        @test isidentical(U2, U)
    end

    # difference
    @test_throws DimensionMismatch difference(U, U3)
    for E2 in (difference(U, U), difference(B, U))
        @test E2 isa EmptySet{N} && dim(E2) == 2
    end
    X = difference(U, B)
    @test X isa UnionSetArray{N,<:HalfSpace} && length(array(X)) == 4 && X == complement(B)

    # distance (between two sets)
    @test_throws DimensionMismatch distance(U, U3)
    @test_throws ArgumentError distance(U, U; p=N(1 // 2))
    for v in (distance(U, U), distance(U, B), distance(B, U), distance(U, Z), distance(Z, U))
        @test v isa N && v == N(0)
    end
    for v in (distance(U, Pe), distance(Pe, U))
        @test v isa N && v == N(Inf)
    end

    # exact_sum
    @test_throws DimensionMismatch exact_sum(U, U3)
    for U2 in (exact_sum(U, U), exact_sum(U, B), exact_sum(B, U))
        @test isidentical(U2, U)
    end

    # intersection
    @test_throws DimensionMismatch intersection(U, U3)
    U2 = intersection(U, U)
    @test isidentical(U2, U)
    for X in (intersection(U, B), intersection(B, U))
        @test X isa BallInf{N} && X == B
    end
    for X in (intersection(U, Pnc), intersection(Pnc, U))
        @test X isa Polygon{N} && X == Pnc
    end

    # isapprox
    @test U ≈ Universe{N}(2)
    @test !(U ≈ U3) && !(U3 ≈ U) && !(U ≈ B) && !(B ≈ U)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(U, U3)
    @test !isdisjoint(U, U)
    res, w = isdisjoint(U, U, true)
    @test !res && w isa Vector{N} && w ∈ U
    @test isdisjoint(U, E) && isdisjoint(E, U) && isdisjoint(U, Pe) && isdisjoint(Pe, U)
    for (res, w) in (isdisjoint(U, E, true), isdisjoint(E, U, true),
                     isdisjoint(U, Pe, true), isdisjoint(Pe, U, true))
        @test res && w isa Vector{N} && isempty(w)
    end
    @test !isdisjoint(U, B) && !isdisjoint(B, U) && !isdisjoint(U, Pnc) && !isdisjoint(Pnc, U)
    for (res, w) in (isdisjoint(U, B, true), isdisjoint(B, U, true), isdisjoint(U, Pnc, true),
                     isdisjoint(Pnc, U, true))
        @test !res && w isa Vector{N} && w ∈ B && w ∈ U
    end

    # isequal
    @test U == Universe{N}(2)
    @test U != U3 && U3 != U && U != B && B != U

    # isequivalent
    @test_throws DimensionMismatch isequivalent(U, U3)
    @test isequivalent(U, U)
    @test !isequivalent(U, B) && !isequivalent(B, U)

    # isstrictsubset
    @test_throws DimensionMismatch U ⊂ U3
    @test !(U ⊂ U)
    res, w = ⊂(U, U, true)
    @test !res && w isa Vector{N} && isempty(w)
    @test !(U ⊂ B)
    res, w = ⊂(U, B, true)
    @test !res && w isa Vector{N} && isempty(w)
    @test B ⊂ U
    res, w = ⊂(B, U, true)
    @test res && w isa Vector{N} && w ∉ B && w ∈ U

    # issubset
    @test_throws DimensionMismatch U ⊆ U3
    for X in (U, B, Pnc)
        @test X ⊆ U
        res, w = ⊆(X, U, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    for X in (B, Pnc)
        @test U ⊈ X
        res, w = ⊆(U, X, true)
        @test !res && w isa Vector{N} && w ∈ U && w ∉ X
    end
    # TODO test `U ⊆ X` with non-Universe `X` for which `isuniversal(X) == true` (currently n/a)

    # linear_combination
    @test_throws DimensionMismatch linear_combination(U, U3)
    for U2 in (linear_combination(U, U),
               linear_combination(U, Pnc), linear_combination(Pnc, U),
               linear_combination(U, B), linear_combination(B, U))
        @test isidentical(U2, U)
    end
    for E2 in (linear_combination(U, Pe), linear_combination(Pe, U))
        @test E2 isa HPolygon{N} && E2 == Pe
    end

    # minkowski_difference
    @test_throws DimensionMismatch minkowski_difference(U, U3)
    # empty difference
    E2 = minkowski_difference(B, U)
    @test E2 isa EmptySet{N} && dim(E2) == 2
    # Universe
    for U2 in (minkowski_difference(U, U), minkowski_difference(U, B), minkowski_difference(U, Z),
               minkowski_difference(LinearMap(N[1 0; 0 1], U), U))
        @test isidentical(U2, U)
    end

    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(U, U3)
    for U2 in (minkowski_sum(U, U), minkowski_sum(U, B), minkowski_sum(B, U))
        @test isidentical(U2, U)
    end
    for X in (minkowski_sum(U, Pe), minkowski_sum(Pe, U))
        @test X isa HPolygon{N} && X == Pe
    end
    for U2 in (minkowski_sum(U, Z), minkowski_sum(Z, U), minkowski_sum(U, B), minkowski_sum(B, U),
               minkowski_sum(U, Pnc), minkowski_sum(Pnc, U))
        @test U2 isa Universe{N} && U2 == U
    end
end

for N in @tN([Float64, Float32])
    U = Universe{N}(2)

    # rationalize
    U2 = rationalize(U)
    @test U2 isa Universe{Rational{Int}} && dim(U2) == 2
    @test_throws MethodError rationalize(U2)

    # exponential_map
    U2 = exponential_map(ones(N, 2, 2), U)
    @test_broken isidentical(U2, U)  # TODO this should change
end
