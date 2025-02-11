function isidentical(::Universe, ::Universe)
    return false
end

function isidentical(U1::Universe{N}, U2::Universe{N}) where {N}
    return U1.dim == U2.dim
end

for _dummy_ in 1:1  # avoid global variable warnings
    # default Float64 constructor
    U = Universe(2)
    @test U isa Universe{Float64}
    @test U.dim == 2
end

for N in [Float64, Float32, Rational{Int}]
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
    @test_throws AssertionError area(U)
    @test_throws AssertionError area(U3)

    # chebyshev_center_radius
    if test_suite_polyhedra
        @test_throws ErrorException chebyshev_center_radius(U)
    end

    # complement
    E2 = complement(U)
    @test E2 isa EmptySet{N} && dim(E2) == 2

    # concretize
    U2 = concretize(U)
    @test isidentical(U, U2)

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
    @test isidentical(U, U2)

    # copy
    U2 = copy(U)
    @test isidentical(U, U2)

    # delaunay
    if isdefined(@__MODULE__, :MiniQhull)  # TODO throwing an error should work without MiniQhull
        @test_throws ArgumentError delaunay(U)
    end

    # diameter
    @test_throws ArgumentError diameter(U)

    # dim
    @test dim(U) == 2

    # eltype
    @test eltype(U) == N
    @test eltype(typeof(U)) == N

    # extrema
    @test extrema(U) == (N[-Inf, -Inf], N[Inf, Inf])
    @test extrema(U, 1) == (N(-Inf), N(Inf))

    # high
    @test high(U) == N[Inf, Inf]
    @test high(U, 1) == N(Inf)

    # isbounded
    @test !isbounded(U)

    # isboundedtype
    @test !isboundedtype(typeof(U))

    # isconvextype
    @test isconvextype(typeof(U))

    # isempty
    @test !isempty(U)

    # isoperation
    @test !isoperation(U)

    # isoperationtype
    @test !isoperationtype(typeof(U))

    # ispolyhedral
    @test ispolyhedral(U)

    # isuniversal
    res, w = isuniversal(U, true)
    @test isuniversal(U) && res && w isa Vector{N} && isempty(w)

    # low
    @test low(U) == N[-Inf, -Inf]
    @test low(U, 1) == N(-Inf)

    # norm
    @test_throws ArgumentError norm(U)

    # polyhedron
    if test_suite_polyhedra
        P = polyhedron(U)
        @test P isa Polyhedra.DefaultPolyhedron
        if N != Float32
            @test P isa Polyhedra.DefaultPolyhedron{N}
        end
        @test size(P.hrep.A) == (0, 2)
    end

    # radius
    @test_throws ArgumentError radius(U)

    # rand
    @test rand(Universe; N=N) isa Universe{N}
    U2 = rand(Universe; N=N, dim=3)
    @test isidentical(U3, U2)

    # rectify
    @test_broken rectify(U)  # TODO see #3687
    # P = rectify(U)
    # Q = HPolyhedron([HalfSpace(N[-1, 0], N(0)), HalfSpace(N[0, -1], N(0))])
    # @test P isa HPolyhedron{N} && isequivalent(P, Q)

    # reflect
    @test reflect(U) == U

    # singleton_list
    @test_throws ArgumentError singleton_list(U)

    # surface
    @test_throws AssertionError surface(U)

    # tosimplehrep
    C, d = tosimplehrep(U)
    @test C isa Matrix{N} && d isa Vector{N} && size(C) == (0, 2) && isempty(d)

    # triangulate
    if test_suite_polyhedra  # TODO this should work without Polyhedra
        @test_throws AssertionError triangulate(U3)
    end

    # vertices_list
    @test_throws ArgumentError vertices_list(U)

    # vertices
    @test_throws ArgumentError vertices(U)

    # volume
    x = volume(U)
    @test x isa N && x == N(Inf)

    # affine_map
    @test_throws AssertionError affine_map(ones(N, 2, 3), U, N[1, 1])
    @test_throws AssertionError affine_map(ones(N, 2, 2), U, N[1])
    if test_suite_polyhedra  # TODO this should work, even without Polyhedra
        @test_broken affine_map(ones(N, 2, 2), U, N[1, 1])
        # U2 = affine_map(ones(N, 2, 2), U, N[1, 1])
        # @test isidentical(U, U2)
        # U2 = affine_map(ones(N, 3, 2), U, N[1, 1, 3])
        # @test isidentical(U3, U2)
    end

    # in
    @test_throws AssertionError N[0] ∈ U
    @test N[0, 0] ∈ U

    # linear_map
    @test_throws AssertionError linear_map(ones(N, 2, 3), U)
    @test_broken linear_map(ones(N, 2, 2), U)  # TODO this should work, even without Polyhedra
    # U2 = linear_map(ones(N, 2, 2), U)
    # @test_broken isidentical(U, U2)
    if test_suite_polyhedra
        @test_broken linear_map(ones(N, 3, 2), U)
        # U2 = linear_map(ones(N, 3, 2), U)
        # @test U2 isa HPolyhedron{N}  # TODO this should change
        # @test_broken isidentical(U3, U2)
    end

    # linear_map_inverse
    U2 = LazySets.linear_map_inverse(ones(N, 2, 3), U)
    @test isidentical(U3, U2)

    # permute
    @test_throws AssertionError permute(U, [1, -1])
    @test_throws AssertionError permute(U, [1, 2, 2])
    U2 = permute(U, [2, 1])
    @test isidentical(U, U2)

    # project
    @test_throws AssertionError project(U, [1, -1])
    @test_throws AssertionError project(U, [1, 2, 3])
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
    @test isidentical(U, U2)
    Z = scale(N(0), U)
    @test Z isa ZeroSet{N} && Z == ZeroSet{N}(2)
    # scale!
    U2 = copy(U)
    scale!(N(2), U2)
    @test isidentical(U, U2)
    @test_throws ArgumentError scale!(N(0), U2)

    # support_function
    @test_throws AssertionError ρ(N[1], U)
    v = ρ(N[-1, 2], U)
    @test v isa N && v == N(Inf)
    v = ρ(N[2, 0], U)
    @test v isa N && v == N(Inf)
    @test ρ(N[0, 0], U) == zero(N)

    # support_vector
    @test_throws AssertionError σ(N[1], U)
    x = σ(N[-1, 2], U)
    @test x isa Vector{N} && x == N[-Inf, Inf]
    x = σ(N[2, 0], U)
    @test x isa Vector{N} && x == N[Inf, 0]
    x = σ(N[0, 0], U)
    @test x isa Vector{N} && x == N[0, 0]

    # translate
    @test_throws AssertionError translate(U, N[1])
    U2 = translate(U, N[1, 2])
    @test isidentical(U, U2)
    # translate!
    @test_throws AssertionError translate!(U, N[1])
    U2 = copy(U)
    translate!(U2, N[1, 2])
    @test isidentical(U, U2)

    # cartesian_product
    for U2 in (cartesian_product(U, U3), cartesian_product(U3, U))
        @test U2 isa Universe{N} && dim(U2) == 5
    end

    # convex_hull (binary)
    @test_throws AssertionError convex_hull(U, U3)
    U2 = convex_hull(U, U)
    @test isidentical(U, U2)
    for U2 in (convex_hull(U, Pnc), convex_hull(Pnc, U), convex_hull(U, E), convex_hull(E, U))
        @test isidentical(U, U2)
    end

    # difference
    @test_throws AssertionError difference(B, U3)
    @test_throws AssertionError difference(U3, B)
    for E2 in (difference(U, U), difference(B, U))
        @test E2 isa EmptySet{N} && dim(E2) == 2
    end
    X = difference(U, B)
    @test X isa UnionSetArray{N,<:HalfSpace} && length(array(X)) == 4 && X == complement(B)

    # distance (between sets)
    @test_throws AssertionError distance(U, U3)
    @test_throws AssertionError distance(U3, U)
    for v in (distance(U, U), distance(U, B), distance(B, U), distance(U, Z), distance(Z, U))
        @test v isa N && v == N(0)
    end
    for v in (distance(U, Pe), distance(Pe, U))
        @test v isa N && v == N(Inf)
    end

    # exact_sum
    @test_throws AssertionError exact_sum(U, U3)
    @test_throws AssertionError exact_sum(U3, U)
    for U2 in (exact_sum(U, U), exact_sum(U, B), exact_sum(B, U))
        @test isidentical(U, U2)
    end

    # intersection
    @test_throws AssertionError intersection(U, U3)
    U2 = intersection(U, U)
    @test isidentical(U, U2)
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
    @test_throws AssertionError isdisjoint(U, U3)
    @test !isdisjoint(U, U)
    res, w = isdisjoint(U, U, true)
    @test !res && w isa Vector{N} && w ∈ U
    @test isdisjoint(U, E) && isdisjoint(E, U) && isdisjoint(U, Pe) && isdisjoint(Pe, U)
    for (res, w) in (isdisjoint(U, E, true), isdisjoint(E, U, true),
                     isdisjoint(U, Pe, true), isdisjoint(Pe, U, true))
        @test res && w isa Vector{N} && w == N[]
    end
    @test !isdisjoint(U, B) && !isdisjoint(B, U) && !isdisjoint(U, Pnc) && !isdisjoint(Pnc, U)
    # TODO add `isdisjoint(U, Pnc, true), isdisjoint(Pnc, U, true)` below once witness production is supported by Polygon
    for (res, w) in (isdisjoint(U, B, true), isdisjoint(B, U, true))
        @test !res && w isa Vector{N} && w ∈ B && w ∈ U
    end

    # isequal
    @test U == Universe{N}(2)
    @test U != U3 && U3 != U && U != B && B != U

    # isequivalent
    @test_throws AssertionError isequivalent(U, U3)
    @test_throws AssertionError isequivalent(U3, U)
    @test isequivalent(U, U)
    @test !isequivalent(U, B) && !isequivalent(B, U)

    # isstrictsubset
    @test_throws AssertionError B ⊂ U3
    @test_throws AssertionError U3 ⊂ B
    @test !(U ⊂ U)
    res, w = ⊂(U, U, true)
    @test !res && w isa Vector{N} && w == N[]
    @test !(U ⊂ B)
    res, w = ⊂(U, B, true)
    @test !res && w isa Vector{N} && w == N[]
    @test B ⊂ U
    res, w = ⊂(B, U, true)
    @test res && w isa Vector{N} && w ∉ B && w ∈ U

    # issubset
    @test_throws AssertionError B ⊆ U3
    @test_throws AssertionError U3 ⊆ B
    for X in (U, B, Pnc)
        @test X ⊆ U
        res, w = ⊆(X, U, true)
        @test res && w isa Vector{N} && w == N[]
    end
    for X in (B, Pnc)
        @test U ⊈ X
        if X === B  # TODO remove once witness production is supported for Polygon
            res, w = ⊆(U, X, true)
            @test !res && w isa Vector{N} && w ∉ X && w ∈ U
        end
    end
    # TODO test `U ⊆ X` with non-Universe `X` for which `isuniversal(X) == true` (currently n/a)

    # linear_combination
    @test_throws AssertionError linear_combination(U, U3)
    for U2 in (linear_combination(U, Pnc), linear_combination(Pnc, U),
               linear_combination(U, B), linear_combination(B, U))
        @test isidentical(U, U2)
    end
    for E2 in (linear_combination(U, Pe), linear_combination(Pe, U))
        @test E2 isa HPolygon{N} && E2 == Pe
    end

    # minkowski_difference
    @test_throws AssertionError minkowski_difference(B, U3)
    @test_throws AssertionError minkowski_difference(U3, B)
    for U2 in (minkowski_difference(U, U), minkowski_difference(U, B), minkowski_difference(U, Z))
        @test isidentical(U, U2)
    end
    E2 = minkowski_difference(B, U)
    @test E2 isa EmptySet{N} && dim(E2) == 2
    # TODO test `minkowski_difference(X, U)` with non-Universe `X` for which `isuniversal(X) == true` (currently n/a)

    # minkowski_sum
    @test_throws AssertionError minkowski_sum(U, U3)
    @test_throws AssertionError minkowski_sum(U3, U)
    for U2 in (minkowski_sum(U, U), minkowski_sum(U, B), minkowski_sum(B, U))
        @test isidentical(U, U2)
    end
    for X in (minkowski_sum(U, Pe), minkowski_sum(Pe, U))
        @test X isa HPolygon{N} && X == Pe
    end
    for U2 in (minkowski_sum(U, Z), minkowski_sum(Z, U), minkowski_sum(U, B), minkowski_sum(B, U),
               minkowski_sum(U, Pnc), minkowski_sum(Pnc, U))
        @test U2 isa Universe{N} && U2 == U
    end
end

for N in [Float64, Float32]
    U = Universe{N}(2)

    # rationalize
    U2 = rationalize(U)
    @test U2 isa Universe{Rational{Int}} && dim(U2) == 2
    @test_throws MethodError rationalize(U2)

    # exponential_map
    @test_throws AssertionError exponential_map(ones(N, 2, 3), U)
    @test_throws AssertionError exponential_map(ones(N, 3, 2), U)
    U2 = exponential_map(ones(N, 2, 2), U)
    @test_broken isidentical(U, U2)  # TODO this should change

    # is_interior_point
    @test_throws AssertionError is_interior_point(N[0], U)
    @test is_interior_point(N[0, 0], U)
end
