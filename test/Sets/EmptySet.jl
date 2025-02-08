function isidentical(::EmptySet{N1}, ::EmptySet{N2}) where {N1,N2}
    return false
end

function isidentical(E1::EmptySet{N}, E2::EmptySet{N}) where {N}
    return E1.dim == E2.dim
end

for _dummy_ in 1:1  # avoid global variable warnings
    # default Float64 constructor
    for E in (EmptySet(2), ∅(2))
        @test E isa EmptySet{Float64}
        @test E.dim == 2
    end
end

for N in [Float64, Float32, Rational{Int}]
    # auxiliary sets
    B = BallInf(ones(N, 2), N(1))
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
    res = area(E)
    @test res isa N && res == N(0)
    @test_throws AssertionError area(E3)

    # complement
    U = complement(E)
    @test U isa Universe{N} && dim(U) == 2

    # concretize
    E2 = concretize(E)
    @test isidentical(E, E2)

    # constrained_dimensions
    @test_broken constrained_dimensions(E)  # TODO implement

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

    # delaunay
    if isdefined(@__MODULE__, :MiniQhull)  # TODO this should throw a normal error without MiniQhull
        @test_broken delaunay(E)
    end

    # diameter
    @test_throws ArgumentError diameter(E)  # TODO this should maybe change
    # res = diameter(E)
    # @test res isa N && res == N(0)

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

    # isoperation
    @test !isoperation(E)

    # isoperationtype
    @test !isoperationtype(typeof(E))

    # ispolyhedral
    @test !ispolyhedral(E)  # TODO this should maybe change

    # isuniversal
    res, w = isuniversal(E, true)
    @test !isuniversal(E) && !res && w isa Vector{N} && w ∉ E

    # low
    @test_throws ArgumentError low(E)
    @test_throws ArgumentError low(E, 1)

    # norm
    @test_throws ArgumentError norm(E)  # TODO this should maybe change
    # res = norm(E)
    # @test res isa N && res == N(0)

    # polyhedron
    if test_suite_polyhedra
        @test_throws MethodError polyhedron(E)  # TODO this should maybe change
    end

    # radius
    @test_throws ArgumentError radius(E)  # TODO this should maybe change
    # res = radius(E)
    # @test res isa N && res == N(0)

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

    # surface
    res = surface(E)
    @test res isa N && res == N(0)

    # tosimplehrep
    @test_throws MethodError tosimplehrep(E)  # TODO this should maybe change

    # triangulate
    if test_suite_polyhedra
        @test_throws AssertionError triangulate(E3)  # TODO this should maybe change
    end

    # vertices_list
    vs = vertices_list(E)
    @test vs == Vector{Vector{N}}() && typeof(vs) == Vector{Vector{N}}

    # vertices
    vs = collect(vertices(E))
    @test vs == Vector{Vector{N}}() && typeof(vs) == Vector{Vector{N}}

    # volume
    @test volume(E) == zero(N)

    # affine_map
    @test_throws AssertionError affine_map(ones(N, 2, 3), E, N[1, 1])
    @test_throws AssertionError affine_map(ones(N, 2, 2), E, N[1])
    E2 = affine_map(ones(N, 2, 2), E, N[1, 1])
    @test isidentical(E, E2)
    E2 = affine_map(ones(N, 3, 2), E, N[1, 1, 3])
    @test isidentical(E3, E2)

    # distance (between point and set)
    @test_throws AssertionError distance(E, N[0])
    x = N[0, 0]
    for res in (distance(E, x), distance(x, E))
        @test res isa N && res == N(Inf)
    end

    # exponential_map / linear_map
    for f in (exponential_map, linear_map)
        @test_throws AssertionError f(ones(N, 2, 3), E)
        E2 = f(ones(N, 2, 2), E)
        @test isidentical(E, E2)
    end
    @test_broken exponential_map(ones(N, 3, 2), E) isa ArgumentError  # TODO rectangular matrix should not be accepted
    E2 = linear_map(ones(N, 3, 2), E)
    @test isidentical(E3, E2)

    # in
    @test_throws AssertionError N[0] ∈ E
    @test N[0, 0] ∉ E

    # linear_map_inverse
    @test_broken LazySets.linear_map_inverse(ones(N, 2, 3), E)  # TODO this should maybe change
    # E2 = LazySets.linear_map_inverse(ones(N, 2, 3), E)
    # @test isidentical(E3, E2)

    # permute
    @test_throws AssertionError permute(E, [1, -1])
    @test_throws AssertionError permute(E, [1, 2, 2])
    E2 = permute(E, [2, 1])
    @test isidentical(E, E2)

    # project
    @test_throws AssertionError project(E, [1, -1])
    @test_throws AssertionError project(E, [1, 2, 3])
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
    @test_throws AssertionError ρ(N[1], E)
    for v in (N[-1, 2], N[2, 0], N[0, 0])
        @test_throws ArgumentError ρ(v, E)
    end

    # support_vector
    @test_throws AssertionError σ(N[1], E)
    for v in (N[-1, 2], N[2, 0], N[0, 0])
        @test_throws ArgumentError σ(v, E)
    end

    # translate
    @test_throws AssertionError translate(E, N[1])
    E2 = translate(E, N[1, 2])
    @test isidentical(E, E2)
    # translate!
    @test_throws AssertionError translate!(E, N[1])
    E2 = copy(E)
    translate!(E2, N[1, 2])
    @test isidentical(E, E2)

    # cartesian_product
    for E2 in (cartesian_product(E, E3), cartesian_product(E3, E))
        @test E2 isa EmptySet{N} && dim(E2) == 5
    end

    # convex_hull (binary)
    @test_throws AssertionError convex_hull(E, E3)
    E2 = convex_hull(E, E)
    @test isidentical(E, E2)
    for X in (convex_hull(E, Pnc), convex_hull(Pnc, E))
        @test X isa LazySet{N} && isequivalent(X, Pc)
    end

    # difference
    @test_throws AssertionError difference(B, E3)
    @test_throws AssertionError difference(E3, B)
    for E2 in (difference(E, E), difference(E, B))
        @test isidentical(E, E2)
    end
    X = difference(B, E)
    @test X isa BallInf{N} && X == B

    # distance (between sets)
    @test_throws AssertionError distance(E, E3)
    @test_throws AssertionError distance(E3, E)
    res = distance(E, E)
    @test res isa N && res == N(Inf)

    # exact_sum / minkowski_sum
    for f in (exact_sum, minkowski_sum)
        @test_throws AssertionError f(E, E3)
        @test_throws AssertionError f(E3, E)
        for E2 in (f(E, E), f(E, B), f(B, E))
            @test isidentical(E, E2)
        end
    end

    # intersection
    @test_throws AssertionError intersection(E, E3)
    for E2 in (intersection(E, E), intersection(E, B), intersection(B, E))
        @test isidentical(E, E2)
    end

    # isapprox
    @test E ≈ E
    @test !(E ≈ E3)
    @test !(E3 ≈ E)

    # isdisjoint
    @test_throws AssertionError isdisjoint(E, E3)
    @test isdisjoint(E, B) && isdisjoint(B, E) && isdisjoint(E, E)
    for (res, w) in (isdisjoint(E, B, true), isdisjoint(B, E, true), isdisjoint(E, E, true))
        @test res && w isa Vector{N} && w == N[]
    end

    # isequal
    @test E == E
    @test !(E == E3) && !(E3 == E)

    # isequivalent
    @test isequivalent(E, E)
    @test_throws AssertionError isequivalent(E, E3)
    @test_throws AssertionError isequivalent(E3, E)

    # isstrictsubset
    @test_throws AssertionError B ⊂ E3
    @test_throws AssertionError E3 ⊂ B
    for X in (E, B)
        @test !(X ⊂ E)
        res, w = ⊂(X, E, true)
        @test !res && w isa Vector{N} && w == N[]
    end
    @test E ⊂ B
    res, w = ⊂(E, B, true)
    @test res && w isa Vector{N} && w ∉ E && w ∈ B

    # issubset
    @test_throws AssertionError B ⊆ E3
    @test_throws AssertionError E3 ⊆ B
    for X in (E, B)
        @test E ⊆ X
        res, w = ⊆(E, X, true)
        @test res && w isa Vector{N} && w == N[]
    end
    @test B ⊈ E
    res, w = ⊆(B, E, true)
    @test !res && w isa Vector{N} && w ∉ E && w ∈ B
    @test Pe ⊆ E
    res, w = ⊆(Pe, E, true)
    @test res && w isa Vector{N} && w == N[]

    # linear_combination
    @test_throws AssertionError linear_combination(E, E3)
    for E2 in (linear_combination(E, Pnc), linear_combination(Pnc, E))
        @test isidentical(E, E2)
    end

    # minkowski_difference
    @test_throws AssertionError minkowski_difference(B, E3)
    @test_throws AssertionError minkowski_difference(E3, B)
    for E2 in (minkowski_difference(E, E), minkowski_difference(E, B))
        @test isidentical(E, E2)
    end
    X = minkowski_difference(B, E)
    @test X isa BallInf{N} && X == B

    # chebyshev_center_radius
    @test_throws ArgumentError chebyshev_center_radius(E)
end

for N in [Float64, Float32]
    E = EmptySet{N}(2)

    # rationalize
    E2 = rationalize(E)
    @test E2 isa EmptySet{Rational{Int}} && dim(E2) == 2
    @test_throws MethodError rationalize(E2)

    # is_interior_point
    @test_throws AssertionError is_interior_point(N[0], E)
    @test !is_interior_point(N[0, 0], E)
end
