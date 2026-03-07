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
    U = @inferred Universe{N}(2)
    @test U isa Universe{N}
    @test U.dim == 2
    U3 = Universe{N}(3)
    @test U3 isa Universe{N}
    @test U3.dim == 3

    # an_element
    x = @inferred an_element(U)
    @test x isa Vector{N} && length(x) == 2

    # area
    @test_throws DimensionMismatch area(Universe{N}(1))
    for X in (U, U3)
        res = @inferred area(X)
        @test res isa N && res == N(Inf)
    end

    # chebyshev_center_radius
    @test_throws ArgumentError chebyshev_center_radius(U)

    # complement
    E2 = @inferred complement(U)
    @test E2 isa EmptySet{N} && dim(E2) == 2

    # concretize
    U2 = @inferred concretize(U)
    @test isidentical(U2, U)

    # constrained_dimensions
    x = @inferred constrained_dimensions(U)
    @test x isa Vector{Int} && isempty(x)

    # constraints_list
    clist = @inferred constraints_list(U)
    @test isempty(clist) && eltype(clist) <: HalfSpace{N}

    # constraints
    citer = @inferred constraints(U)
    @test isempty(citer) && eltype(citer) <: AbstractVector{N}

    # convex_hull (unary)
    U2 = @inferred convex_hull(U)
    @test isidentical(U2, U)

    # copy
    U2 = @inferred copy(U)
    @test isidentical(U2, U)

    # diameter
    @test_throws ArgumentError diameter(U, N(1 // 2))
    @test_throws ArgumentError diameter(U)
    @test_throws ArgumentError diameter(U, Inf)
    @test_throws ArgumentError diameter(U, 2)

    # dim
    @test (@inferred dim(U)) == 2
    @test (@inferred dim(U3)) == 3

    # eltype
    @test (@inferred eltype(U)) == N
    @test (@inferred eltype(typeof(U))) == N

    # extrema
    @test (@inferred extrema(U)) == (N[-Inf, -Inf], N[Inf, Inf])
    @test_throws DimensionMismatch extrema(U, 3)
    @test (@inferred extrema(U, 1)) == (N(-Inf), N(Inf))

    # high
    @test (@inferred high(U)) == N[Inf, Inf]
    @test_throws DimensionMismatch high(U, 3)
    @test (@inferred high(U, 1)) == N(Inf)

    # isbounded
    @test !(@inferred isbounded(U))

    # isboundedtype
    @test !(@inferred isboundedtype(typeof(U)))

    # isconvex
    @test @inferred isconvex(U)

    # isconvextype
    @test @inferred isconvextype(typeof(U))

    # isempty
    @test !(@inferred isempty(U))
    @test_broken @inferred isempty(U, true)  # TODO make this type-stable
    res, w = isempty(U, true)
    @test !res && w isa Vector{N} && w ∈ U

    # isoperation
    @test !(@inferred isoperation(U))

    # isoperationtype
    @test !(@inferred isoperationtype(typeof(U)))

    # ispolyhedral
    @test @inferred ispolyhedral(U)

    # ispolyhedraltype
    @test @inferred ispolyhedraltype(typeof(U))

    # ispolytopic
    @test !(@inferred ispolytopic(U))

    # ispolytopictype
    @test !(@inferred ispolytopictype(typeof(U)))

    # isuniversal
    @test @inferred isuniversal(U)
    @test_broken @inferred isuniversal(U, true)  # TODO make this type-stable
    res, w = isuniversal(U, true)
    @test res && w isa Vector{N} && isempty(w)

    # low
    @test (@inferred low(U)) == N[-Inf, -Inf]
    @test_throws DimensionMismatch low(U, 3)
    @test (@inferred low(U, 1)) == N(-Inf)

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
    U2 = rand(Universe; N=N)
    @test isidentical(U2, U)
    U2 = rand(Universe; N=N, dim=3)
    @test isidentical(U2, U3)

    # rectify
    P = @inferred rectify(U)
    Q = HPolyhedron([HalfSpace(N[-1, 0], N(0)), HalfSpace(N[0, -1], N(0))])
    @test P isa HPolyhedron{N} && isequivalent(P, Q)

    # reflect
    U2 = @inferred reflect(U)
    @test U2 == U

    # singleton_list
    @test_throws ArgumentError singleton_list(U)

    # tosimplehrep
    C, d = @inferred tosimplehrep(U)
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
    x = @inferred volume(U)
    @test x isa N && x == N(Inf)

    # affine_map (part 1)
    @test_throws DimensionMismatch affine_map(ones(N, 2, 3), U, N[1, 1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 2), U, N[1])
    U2 = affine_map(N[1 0; 0 1], U, N[1, 1])
    @test isidentical(U2, U)

    # distance (between point and set)
    @test_throws DimensionMismatch distance(U, N[0])
    @test_throws ArgumentError distance(U, N[0, 0]; p=N(1 // 2))
    x = N[0, 0]
    for res in ((@inferred distance(U, x)), @inferred distance(x, U))
        @test res isa N && res == N(0)
    end

    # exponential_map
    @test_throws DimensionMismatch exponential_map(ones(N, 1, 1), U)
    @test_throws DimensionMismatch exponential_map(ones(N, 3, 2), U)

    # in
    @test_throws DimensionMismatch N[0] ∈ U
    @test @inferred N[0, 0] ∈ U

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0], U)
    @test_throws ArgumentError is_interior_point(N[0, 0], U; ε=N(0))
    @test_throws ArgumentError is_interior_point(N[0, 0], U; p=N(1 // 2))
    if N <: AbstractFloat
        @static if VERSION >= v"1.11"
            @test @inferred is_interior_point(N[0, 0], U)
        else
            @test is_interior_point(N[0, 0], U)
        end
    else
        @test_throws ArgumentError is_interior_point(N[0, 0], U)
        @test @inferred is_interior_point(N[0, 0], U; ε=1 // 100)
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([0.0, 0.0], U)
    end

    # linear_map (part 1)
    @test_throws DimensionMismatch linear_map(ones(N, 2, 1), U)
    U2 = linear_map(N[1 0; 0 1], U)
    @test isidentical(U2, U)

    # linear_map_inverse
    U2 = @inferred LazySets.linear_map_inverse(ones(N, 2, 3), U)
    @test isidentical(U2, U3)

    # permute
    @test_throws DimensionMismatch permute(U, [1])
    @test_throws DimensionMismatch permute(U, [1, -1])
    @test_throws DimensionMismatch permute(U, [1, 3])
    @test_throws ArgumentError permute(U, [1, 1])
    U2 = @inferred permute(U, [1, 2])
    @test isidentical(U2, U)
    U2 = @inferred permute(U, [2, 1])
    @test isidentical(U2, U)

    # project
    @test_throws DimensionMismatch project(U, [1, 2, 3])
    @test_throws DimensionMismatch project(U, [1, -1])
    @test_throws DimensionMismatch project(U, [1, 3])
    @test_throws ArgumentError project(U, [1, 1])
    U2 = @inferred project(U, [2])
    @test U2 isa Universe{N} && dim(U2) == 1

    # sample
    x = @inferred sample(U)
    @test x isa Vector{N} && length(x) == 2
    xs = @inferred sample(U, 2)
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
    @inferred scale!(N(2), U2)
    @test isidentical(U2, U)
    @test_throws ArgumentError scale!(N(0), U2)

    # support_function
    @test_throws DimensionMismatch ρ(N[1], U)
    v = @inferred ρ(N[-1, 2], U)
    @test v isa N && v == N(Inf)
    v = @inferred ρ(N[2, 0], U)
    @test v isa N && v == N(Inf)
    v = @inferred ρ(N[0, 0], U)
    @test v isa N && v == N(0)

    # support_vector
    @test_throws DimensionMismatch σ(N[1], U)
    x = @inferred σ(N[-1, 2], U)
    @test x isa Vector{N} && x == N[-Inf, Inf]
    x = @inferred σ(N[2, 0], U)
    @test x isa Vector{N} && x == N[Inf, 0]
    x = @inferred σ(N[0, 0], U)
    @test x isa Vector{N} && x == N[0, 0]

    # translate
    @test_throws DimensionMismatch translate(U, N[1])
    U2 = @inferred translate(U, N[1, 2])
    @test isidentical(U2, U)
    # translate!
    @test_throws DimensionMismatch translate!(U, N[1])
    U2 = copy(U)
    @inferred translate!(U2, N[1, 2])
    @test isidentical(U2, U)

    # cartesian_product
    for U2 in ((@inferred cartesian_product(U, U3)), @inferred cartesian_product(U3, U))
        @test U2 isa Universe{N} && dim(U2) == 5
    end

    # convex_hull (binary)
    @test_throws DimensionMismatch convex_hull(U, U3)
    U2 = @inferred convex_hull(U, U)
    @test isidentical(U2, U)
    for X in (E, Pnc)
        U2 = @inferred convex_hull(U, X)
        @test isidentical(U2, U)
        U2 = @inferred convex_hull(X, U)
        @test isidentical(U2, U)
    end

    # difference
    @test_throws DimensionMismatch difference(U, U3)
    for E2 in ((@inferred difference(U, U)), @inferred difference(B, U))
        @test E2 isa EmptySet{N} && dim(E2) == 2
    end
    X = @inferred difference(U, B)
    @test X isa UnionSetArray{N,<:HalfSpace} && length(array(X)) == 4 && X == complement(B)

    # distance (between two sets)
    @test_throws DimensionMismatch distance(U, U3)
    @test_throws ArgumentError distance(U, U; p=N(1 // 2))
    v = @inferred distance(U, U)
    @test v isa N && v == N(0)
    for X in (B, Z)
        v = @inferred distance(U, X)
        @test v isa N && v == N(0)
        v = @inferred distance(X, U)
        @test v isa N && v == N(0)
    end
    for v in ((@inferred distance(U, Pe)), @inferred distance(Pe, U))
        @test v isa N && v == N(Inf)
    end

    # exact_sum
    @test_throws DimensionMismatch exact_sum(U, U3)
    for U2 in ((@inferred exact_sum(U, U)), (@inferred exact_sum(U, B)), @inferred exact_sum(B, U))
        @test isidentical(U2, U)
    end

    # intersection
    @test_throws DimensionMismatch intersection(U, U3)
    U2 = @inferred intersection(U, U)
    @test isidentical(U2, U)
    for X in ((@inferred intersection(U, B)), @inferred intersection(B, U))
        @test X isa BallInf{N} && X == B
    end
    for X in ((@inferred intersection(U, Pnc)), @inferred intersection(Pnc, U))
        @test X isa Polygon{N} && X == Pnc
    end

    # isapprox
    @test @inferred U ≈ Universe{N}(2)
    @test !(@inferred U ≈ U3) && !(@inferred U3 ≈ U) && !(@inferred U ≈ B) && !(@inferred B ≈ U)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(U, U3)
    @static if VERSION >= v"1.12"
        @test_broken !(@inferred isdisjoint(U, U))  # TODO make this type-stable
        @test !isdisjoint(U, U)
    else
        @test !(@inferred isdisjoint(U, U))
    end
    @static if VERSION >= v"1.11"
        @test_broken @inferred isdisjoint(U, U, true)  # TODO make this type-stable
    end
    res, w = isdisjoint(U, U, true)
    @test !res && w isa Vector{N} && w ∈ U
    @test (@inferred isdisjoint(U, E)) && (@inferred isdisjoint(E, U))
    @static if VERSION >= v"1.12"
        @test_broken @inferred isdisjoint(U, Pe)  # TODO make this type-stable
    end
    @test isdisjoint(U, Pe) && isdisjoint(Pe, U)
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
    @test @inferred U == Universe{N}(2)
    @test (@inferred U != U3) && (@inferred U3 != U) && (@inferred U != B) && @inferred B != U

    # isequivalent
    @test_throws DimensionMismatch isequivalent(U, U3)
    @test @inferred isequivalent(U, U)
    @test_broken @inferred isequivalent(B, U)  # TODO make this type-stable
    @test !(@inferred isequivalent(U, B)) && !isequivalent(B, U)

    # isstrictsubset
    @test_throws DimensionMismatch U ⊂ U3
    @test_broken !(@inferred U ⊂ U)  # TODO make this type-stable
    @test !(U ⊂ U)
    @test_broken @inferred ⊂(U, U, true)  # TODO make this type-stable
    res, w = ⊂(U, U, true)
    @test !res && w isa Vector{N} && isempty(w)
    @test !(U ⊂ B)
    res, w = ⊂(U, B, true)
    @test !res && w isa Vector{N} && w ∈ U && w ∉ B
    @test_broken @inferred B ⊂ U  # TODO make this type-stable
    @test B ⊂ U
    res, w = ⊂(B, U, true)
    @test res && w isa Vector{N} && w ∉ B && w ∈ U

    # issubset
    @test_throws DimensionMismatch U ⊆ U3
    @test_broken @inferred ⊆(U, U, true)  # TODO make this type-stable
    for X in (U, B, Pnc)
        @test @inferred X ⊆ U
        res, w = ⊆(X, U, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    for X in (B, Pnc)
        @test @inferred U ⊈ X
        res, w = ⊆(U, X, true)
        @test !res && w isa Vector{N} && w ∈ U && w ∉ X
    end
    U2 = Universe{N}(1)
    X = Line(N[0], N[1])
    @static if VERSION >= v"1.12"
        @test_broken @inferred U2 ⊆ X  # TODO make this type-stable
    end
    @test U2 ⊆ X
    res, w = ⊆(U2, X, true)
    @test res && w isa Vector{N} && isempty(w)

    # linear_combination
    @test_throws DimensionMismatch linear_combination(U, U3)
    for U2 in ((@inferred linear_combination(U, U)), (@inferred linear_combination(U, B)),
               @inferred linear_combination(B, U))
        @test isidentical(U2, U)
    end
    for U2 in (linear_combination(U, Pnc), linear_combination(Pnc, U))
        @test isidentical(U2, U)
    end
    for E2 in (linear_combination(U, Pe), linear_combination(Pe, U))
        @test E2 isa HPolygon{N} && E2 == Pe
    end

    # minkowski_difference
    @test_throws DimensionMismatch minkowski_difference(U, U3)
    # empty difference
    E2 = @inferred minkowski_difference(B, U)
    @test E2 isa EmptySet{N} && dim(E2) == 2
    # Universe
    for X in (U, B, Z)
        U2 = @inferred minkowski_difference(U, X)
        @test isidentical(U2, U)
    end
    X = LinearMap(N[1 0; 0 1], U)
    U2 = minkowski_difference(X, U)
    @test isidentical(U2, U)

    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(U, U3)
    for U2 in ((@inferred minkowski_sum(U, U)), (@inferred minkowski_sum(U, B)),
               @inferred minkowski_sum(B, U))
        @test isidentical(U2, U)
    end
    for X in (minkowski_sum(U, Pe), minkowski_sum(Pe, U))
        @test X isa HPolygon{N} && X == Pe
    end
    for X in (Z, B)
        for U2 in ((@inferred minkowski_sum(U, X)), @inferred minkowski_sum(X, U))
            @test U2 isa Universe{N} && U2 == U
        end
    end
    for U2 in (minkowski_sum(U, Pnc), minkowski_sum(Pnc, U))
        @test U2 isa Universe{N} && U2 == U
    end
end

for N in @tN([Float64, Float32])
    U = Universe{N}(2)

    # rationalize
    U2 = @inferred rationalize(U)
    @test U2 isa Universe{Rational{Int}} && dim(U2) == 2
    @test_throws MethodError rationalize(U2)

    # affine_map (part 2)
    X = affine_map(N[1 0; 0 1; 0 0], U, N[1, 1, 3])
    @test X isa HPolyhedron{N} && isequivalent(X, Hyperplane(N[0, 0, 1], N(3)))

    # exponential_map
    U2 = @inferred exponential_map(ones(N, 2, 2), U)
    @test isidentical(U2, U)

    # linear_map (part 2)
    X = linear_map(N[1 0; 0 1; 0 0], U)
    @test X isa HPolyhedron{N} && isequivalent(X, Hyperplane(N[0, 0, 1], N(0)))
end

for N in [Float64]
    U = Universe{N}(2)

    # affine_map (part 3)
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :CDDLib)
        X = affine_map(N[1 2; 0 0], U, N[1, 1])  # projection to axis
        @test X isa HPolyhedron{N} && isequivalent(X, Hyperplane(N[0, 1], N(1)))
        X = affine_map(ones(N, 2, 2), U, N[2, 0])  # projection to line
        @test X isa HPolyhedron{N} && isequivalent(X, Line2D(N[1, -1], N(2)))
        X = affine_map(zeros(N, 2, 2), U, N[2, 0])  # zero map
        @test X isa HPolyhedron{N} && isequivalent(X, Singleton(N[2, 0]))
    end

    # linear_map (part 3)
    @static if isdefined(@__MODULE__, :Polyhedra) && isdefined(@__MODULE__, :CDDLib)
        X = linear_map(N[1 2; 0 0], U)  # projection to axis
        @test X isa HPolyhedron{N} && isequivalent(X, Hyperplane(N[0, 1], N(0)))
        X = linear_map(ones(N, 2, 2), U)  # projection to line
        @test X isa HPolyhedron{N} && isequivalent(X, Line2D(N[1, -1], N(0)))
        X = linear_map(zeros(N, 2, 2), U)  # zero map
        @test X isa HPolyhedron{N} && isequivalent(X, ZeroSet{N}(2))
    end
end
