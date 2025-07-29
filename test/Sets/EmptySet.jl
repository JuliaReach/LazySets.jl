using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

function isidentical(::EmptySet, ::EmptySet)
    return false
end

function isidentical(E1::EmptySet{N}, E2::EmptySet{N}) where {N}
    return E1.dim == E2.dim
end

let
    # default Float64 constructor
    for E in (EmptySet(2), ∅(2))
        @test E isa EmptySet{Float64}
        @test E.dim == 2
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # auxiliary sets
    B = BallInf(ones(N, 2), N(1))
    U = Universe{N}(2)
    Z = ZeroSet{N}(2)
    B1 = Ball1(ones(N, 2), N(1))
    Pe = HPolygon([HalfSpace(N[1, 0], N(0)), HalfSpace(N[-1, 0], N(-1)),  # empty
                   HalfSpace(N[0, 1], N(0)), HalfSpace(N[0, -1], N(0))])
    Pnc = Polygon([N[0, 0], N[3, 0], N[1, 1], N[0, 3]])  # nonconvex
    Pc = VPolygon([N[0, 0], N[3, 0], N[0, 3]])  # convex hull of `Pnc`

    # constructor
    E = EmptySet{N}(2)
    @test E isa EmptySet{N}
    @test E.dim == 2
    E3 = EmptySet{N}(3)
    @test E3 isa EmptySet{N}
    @test E3.dim == 3

    # convert
    E2 = convert(EmptySet, Pe)
    @test isidentical(E, E2)
    @test_throws AssertionError convert(EmptySet, B)

    # an_element
    @test_throws ArgumentError an_element(E)

    # area
    @test_throws DimensionMismatch area(EmptySet{N}(1))
    for X in (E, E3)
        res = area(X)
        @test res isa N && res == N(0)
    end

    # chebyshev_center_radius
    @test_throws ArgumentError chebyshev_center_radius(E)

    # complement
    U2 = complement(E)
    @test U2 isa Universe{N} && dim(U2) == 2

    # concretize
    E2 = concretize(E)
    @test isidentical(E, E2)

    # constrained_dimensions
    @test constrained_dimensions(E) == 1:2

    # constraints_list
    @test_throws MethodError constraints_list(E)  # TODO this should maybe change

    # constraints
    @test_throws MethodError constraints(E)  # TODO this should maybe change

    # convex_hull (unary)
    E2 = convex_hull(E)
    @test isidentical(E, E2)

    # copy
    E2 = copy(E)
    @test isidentical(E, E2)

    # diameter
    for res in (diameter(E), diameter(E, Inf), diameter(E, 2))
        @test res isa N && res == N(0)
    end

    # dim
    @test dim(E) == 2

    # eltype
    @test eltype(E) == N
    @test eltype(typeof(E)) == N

    # extrema
    @test_throws ArgumentError extrema(E)
    @test_throws ArgumentError extrema(E, 1)

    # high
    @test_throws ArgumentError high(E)
    @test_throws ArgumentError high(E, 1)

    # isbounded
    @test isbounded(E)

    # isboundedtype
    @test isboundedtype(typeof(E))

    # isconvextype
    @test isconvextype(typeof(E))

    # isempty
    @test isempty(E)
    @test_broken isempty(E, true)  # TODO this should change

    # isoperation
    @test !isoperation(E)

    # isoperationtype
    @test !isoperationtype(typeof(E))

    # ispolyhedral
    @test !ispolyhedral(E)  # TODO this should maybe change

    # isuniversal
    @test !isuniversal(E)
    res, w = isuniversal(E, true)
    @test !res && w isa Vector{N} && w ∉ E

    # low
    @test_throws ArgumentError low(E)
    @test_throws ArgumentError low(E, 1)

    # norm
    for res in (norm(E), norm(E, Inf), norm(E, 2))
        @test res isa N && res == N(0)
    end

    # polyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test_throws MethodError polyhedron(E)  # TODO this should maybe change
    end

    # radius
    for res in (radius(E), radius(E, Inf), radius(E, 2))
        @test res isa N && res == N(0)
    end

    # rand
    E2 = rand(EmptySet; N=N)
    @test isidentical(E, E2)
    E2 = rand(EmptySet; N=N, dim=3)
    @test isidentical(E3, E2)

    # rectify
    @test rectify(E) == E

    # reflect
    @test reflect(E) == E

    # singleton_list
    res = singleton_list(E)
    T = VERSION < v"1.7" ? Singleton : Singleton{N,Vector{N}}
    @test res isa Vector{T} && isempty(res)

    # tosimplehrep
    @test_throws MethodError tosimplehrep(E)  # TODO this should maybe change

    # triangulate
    @test_throws ArgumentError triangulate(E)

    # triangulate_faces
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test_throws AssertionError triangulate_faces(E3)  # TODO this should maybe change
    end

    # vertices_list
    res = vertices_list(E)
    @test res isa Vector{Vector{N}} && isempty(res)

    # vertices
    res = collect(vertices(E))
    @test res isa Vector{Vector{N}} && isempty(res)

    # volume
    @test volume(E) == N(0)

    # affine_map
    @test_throws DimensionMismatch affine_map(ones(N, 2, 3), E, N[1, 1])
    @test_throws DimensionMismatch affine_map(ones(N, 2, 2), E, N[1])
    E2 = affine_map(ones(N, 2, 2), E, N[1, 1])
    @test isidentical(E, E2)
    E2 = affine_map(ones(N, 3, 2), E, N[1, 1, 3])
    @test isidentical(E3, E2)

    # distance (between point and set)
    @test_throws DimensionMismatch distance(E, N[0])
    x = N[0, 0]
    for res in (distance(E, x), distance(x, E))
        @test res isa N && res == N(Inf)
    end

    # exponential_map
    @test_throws DimensionMismatch exponential_map(ones(N, 1, 1), E)
    @test_throws DimensionMismatch exponential_map(ones(N, 3, 2), E)
    E2 = exponential_map(ones(N, 2, 2), E)
    @test isidentical(E, E2)

    # in
    @test_throws DimensionMismatch N[0] ∈ E
    @test N[0, 0] ∉ E

    # is_interior_point
    @test_throws DimensionMismatch is_interior_point(N[0], E)
    @test_throws DimensionMismatch is_interior_point(N[0], E; ε=N(0))
    if N <: AbstractFloat
        @test !is_interior_point(N[0, 0], E)
    else
        @test_throws ArgumentError is_interior_point(N[0, 0], E)
        @test !is_interior_point(N[0, 0], E; ε=1 // 100)
        # incompatible numeric type
        @test_throws ArgumentError is_interior_point([0.0, 0.0], E)
    end

    # linear_map
    @test_throws DimensionMismatch linear_map(ones(N, 2, 1), E)
    E2 = linear_map(ones(N, 2, 2), E)
    @test isidentical(E, E2)
    E2 = linear_map(ones(N, 3, 2), E)
    @test isidentical(E3, E2)

    # linear_map_inverse
    E2 = LazySets.linear_map_inverse(ones(N, 2, 3), E)
    @test isidentical(E3, E2)

    # permute
    @test_throws DimensionMismatch permute(E, [1])
    @test_throws DimensionMismatch permute(E, [1, -1])
    @test_throws DimensionMismatch permute(E, [1, 3])
    E2 = permute(E, [2, 1])
    @test isidentical(E, E2)

    # project
    @test_throws DimensionMismatch project(E, [1, 2, 3])
    @test_throws DimensionMismatch project(E, [1, -1])
    @test_throws DimensionMismatch project(E, [1, 3])
    E2 = project(E, [2])
    @test E2 isa EmptySet{N} && dim(E2) == 1

    # sample
    @test_throws ArgumentError sample(E)
    @test_throws ArgumentError sample(E, 2)

    # scale
    E2 = scale(N(2), E)
    @test isidentical(E, E2)
    E2 = scale(N(0), E)
    @test isidentical(E, E2)
    # scale!
    E2 = copy(E)
    scale!(N(2), E2)
    @test isidentical(E, E2)
    scale!(N(0), E2)
    @test isidentical(E, E2)

    # support_function
    @test_throws DimensionMismatch ρ(N[1], E)
    for x in (N[-1, 2], N[2, 0], N[0, 0])
        @test_throws ArgumentError ρ(x, E)
    end

    # support_vector
    @test_throws DimensionMismatch σ(N[1], E)
    for x in (N[-1, 2], N[2, 0], N[0, 0])
        @test_throws ArgumentError σ(x, E)
    end

    # translate
    @test_throws DimensionMismatch translate(E, N[1])
    E2 = translate(E, N[1, 2])
    @test isidentical(E, E2)
    # translate!
    @test_throws DimensionMismatch translate!(E, N[1])
    E2 = copy(E)
    translate!(E2, N[1, 2])
    @test isidentical(E, E2)

    # cartesian_product
    for E2 in (cartesian_product(E, E3), cartesian_product(E3, E))
        @test E2 isa EmptySet{N} && dim(E2) == 5
    end
    for E2 in (cartesian_product(E, B), cartesian_product(B, E))
        @test E2 isa EmptySet{N} && dim(E2) == 4
    end

    # convex_hull (binary)
    @test_throws DimensionMismatch convex_hull(E, E3)
    E2 = convex_hull(E, E)
    @test isidentical(E, E2)
    for X in (convex_hull(E, Pnc), convex_hull(Pnc, E))
        @test X isa LazySet{N} && isequivalent(X, Pc)
    end

    # difference
    @test_throws DimensionMismatch difference(E, E3)
    @test_throws DimensionMismatch difference(B, E3)
    @test_throws DimensionMismatch difference(E3, B)
    for E2 in (difference(E, E), difference(E, B), difference(E, U))
        @test isidentical(E, E2)
    end
    X = difference(B, E)
    @test X isa BallInf{N} && X == B
    U2 = difference(U, E)
    @test U2 isa Universe{N} && U2 == U

    # distance (between two sets)
    @test_throws DimensionMismatch distance(E, E3)
    for v in (distance(E, E), distance(B1, E), distance(E, B1), distance(U, E), distance(E, U),
              distance(E, B), distance(B, E), distance(E, Z), distance(Z, E))
        @test v isa N && v == N(Inf)
    end

    # exact_sum
    @test_throws DimensionMismatch exact_sum(E, E3)
    for E2 in (exact_sum(E, E), exact_sum(E, B), exact_sum(B, E))
        @test isidentical(E, E2)
    end

    # intersection
    @test_throws DimensionMismatch intersection(E, E3)
    for E2 in (intersection(E, E), intersection(E, B), intersection(B, E),
               intersection(E, U), intersection(U, E), intersection(E, Z), intersection(Z, E))
        @test isidentical(E, E2)
    end

    # isapprox
    @test E ≈ EmptySet{N}(2)
    @test !(E ≈ E3) && !(E3 ≈ E) && !(E ≈ Pe) && !(Pe ≈ E)

    # isdisjoint
    @test_throws DimensionMismatch isdisjoint(E, E3)
    @test isdisjoint(E, E) && isdisjoint(E, B) && isdisjoint(B, E) &&
          isdisjoint(E, Pnc) && isdisjoint(Pnc, E)
    for (res, w) in (isdisjoint(E, E, true), isdisjoint(E, B, true), isdisjoint(B, E, true),
                     isdisjoint(E, Pnc, true), isdisjoint(Pnc, E, true))
        @test res && w isa Vector{N} && isempty(w)
    end

    # isequal
    @test E == EmptySet{N}(2)
    @test E != E3 && E3 != E && E != B && B != E

    # isequivalent
    @test_throws DimensionMismatch isequivalent(E, E3)
    @test isequivalent(E, E)
    @test !isequivalent(E, B) && !isequivalent(B, E)

    # isstrictsubset
    @test_throws DimensionMismatch E ⊂ E3
    @test_throws DimensionMismatch B ⊂ E3
    @test_throws DimensionMismatch E3 ⊂ B
    for X in (E, B)
        @test !(X ⊂ E)
        res, w = ⊂(X, E, true)
        @test !res && w isa Vector{N} && isempty(w)
    end
    @test E ⊂ B
    res, w = ⊂(E, B, true)
    @test res && w isa Vector{N} && w ∉ E && w ∈ B

    # issubset
    @test_throws DimensionMismatch B ⊆ E3
    @test_throws DimensionMismatch E3 ⊆ B
    for X in (E, B, Pnc)
        @test E ⊆ X
        res, w = ⊆(E, X, true)
        @test res && w isa Vector{N} && isempty(w)
    end
    for X in (B, Pnc)
        @test X ⊈ E
        res, w = ⊆(X, E, true)
        @test !res && w isa Vector{N} && w ∈ X && w ∉ E
    end
    @test Pe ⊆ E
    res, w = ⊆(Pe, E, true)
    @test res && w isa Vector{N} && isempty(w)

    # linear_combination
    @test_throws DimensionMismatch linear_combination(E, E3)
    for E2 in (linear_combination(E, Pnc), linear_combination(Pnc, E),
               linear_combination(E, B), linear_combination(B, E),
               linear_combination(E, U), linear_combination(U, E))
        @test isidentical(E, E2)
    end

    # minkowski_difference
    @test_throws DimensionMismatch minkowski_difference(E, E3)
    @test_throws DimensionMismatch minkowski_difference(B, E3)
    @test_throws DimensionMismatch minkowski_difference(E3, B)
    # empty difference
    for E2 in (minkowski_difference(E, E), minkowski_difference(E, B),
               minkowski_difference(E, U), minkowski_difference(E, Z))
        @test isidentical(E, E2)
    end
    # nonempty difference
    X = minkowski_difference(B, E)
    @test X isa BallInf{N} && X == B
    # Universe
    U2 = minkowski_difference(U, E)
    @test U2 isa Universe{N} && dim(U2) == 2

    # minkowski_sum
    @test_throws DimensionMismatch minkowski_sum(E, E3)
    for E2 in (minkowski_sum(E, E), minkowski_sum(E, B), minkowski_sum(B, E), minkowski_sum(U, E),
               minkowski_sum(E, U), minkowski_sum(E, Z), minkowski_sum(Z, E), minkowski_sum(E, B),
               minkowski_sum(B, E))
        @test isidentical(E, E2)
    end
end

for N in @tN([Float64, Float32])
    E = EmptySet{N}(2)

    # rationalize
    E2 = rationalize(E)
    @test E2 isa EmptySet{Rational{Int}} && dim(E2) == 2
    @test_throws MethodError rationalize(E2)
end
