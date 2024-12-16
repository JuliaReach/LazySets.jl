for _dummy_ in 1:1  # avoid global variable warnings
    # default Float64 constructor
    for E in (EmptySet(2), ∅(2))
        @test E isa EmptySet{Float64}
        @test E.dim == 2
    end
end

for N in [Float64, Rational{Int}, Float32]
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

    # convert
    E2 = convert(EmptySet, Pe)
    @test E2 isa EmptySet{N} && dim(E2) == 2
    @test_throws AssertionError convert(EmptySet, B)

    # an_element
    @test_throws ArgumentError an_element(E)

    # area
    res = area(E)
    @test res isa N && res == N(0)

    # complement
    U = complement(E)
    @test U isa Universe{N} && dim(U) == 2

    # concretize
    E2 = concretize(E)
    @test E2 isa EmptySet{N} && dim(E2) == 2

    # constraints_list
    @test_throws MethodError constraints_list(E)  # TODO this should maybe change

    # constraints
    @test_throws MethodError constraints(E)  # TODO this should maybe change

    # copy
    E2 = copy(E)
    @test E2 isa EmptySet{N} && dim(E2) == 2

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
    @test !isuniversal(E) && !res && w isa Vector{N} && length(w) == 2 && w ∉ E

    # low
    @test_throws ArgumentError low(E)
    @test_throws ArgumentError low(E, 1)

    # norm
    @test_throws ArgumentError norm(E)  # TODO this should maybe change
    # res = norm(E)
    # @test res isa N && res == N(0)

    # radius
    @test_throws ArgumentError radius(E)  # TODO this should maybe change
    # res = radius(E)
    # @test res isa N && res == N(0)

    # rand
    @test rand(EmptySet; N=N) isa EmptySet{N}
    E2 = rand(EmptySet; N=N, dim=3)
    @test E2 isa EmptySet{N} && dim(E2) == 3

    # rectify
    @test rectify(E) == E

    # reflect
    @test reflect(E) == E

    # sample
    @test_throws ArgumentError sample(E)
    @test_throws ArgumentError sample(E, 2)

    # singleton_list
    res = singleton_list(E)
    T = VERSION < v"1.7" ? Singleton : Singleton{N, Vector{N}}
    @test res isa Vector{T} && isempty(res)

    # surface
    res = area(E)
    @test res isa N && res == N(0)

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
    @test E2 isa EmptySet{N} && dim(E2) == 2
    E2 = affine_map(ones(N, 3, 2), E, N[1, 1, 3])
    @test E2 isa EmptySet{N} && dim(E2) == 3

    # exponential_map / linear_map
    for f in (exponential_map, linear_map)
        @test_throws AssertionError f(ones(N, 2, 3), E)
        E2 = f(ones(N, 2, 2), E)
        @test E2 isa EmptySet{N} && dim(E2) == 2
        E2 = f(ones(N, 3, 2), E)
        @test E2 isa EmptySet{N} && dim(E2) == 3
    end

    # in
    @test_throws AssertionError N[0] ∈ E
    @test N[0, 0] ∉ E

    # permute
    @test_throws AssertionError permute(E, [1, -1])
    @test_throws AssertionError permute(E, [1, 2, 2])
    E2 = permute(E, [2, 1])
    @test E2 isa EmptySet{N} && dim(E2) == 2

    # project
    @test_throws AssertionError project(E, [1, -1])
    @test_throws AssertionError project(E, [1, 2, 3])
    E2 = project(E, [2])
    @test E2 isa EmptySet{N} && dim(E2) == 1

    # scale
    E2 = scale(N(2), E)
    @test E2 isa EmptySet{N} && dim(E2) == 2
    # scale!
    E3 = EmptySet{N}(2)
    scale!(N(2), E3)
    @test E3 isa EmptySet{N} && dim(E3) == 2

    # support_function
    @test_throws AssertionError ρ(N[1], E)
    @test_throws ArgumentError ρ(N[1, 1], E)

    # support_vector
    @test_throws AssertionError σ(N[1], E)
    @test_throws ArgumentError σ(N[1, 1], E)

    # translate
    @test_throws AssertionError translate(E, N[1])
    E2 = translate(E, N[1, 2])
    @test E2 isa EmptySet{N} && dim(E2) == 2
    # translate!
    @test_throws AssertionError translate!(E, N[1])
    E2 = EmptySet{N}(2)
    translate!(E2, N[1, 2])
    @test E2 isa EmptySet{N} && dim(E2) == 2

    # cartesian_product
    E2 = EmptySet{N}(3)
    for X in (cartesian_product(E, E2), cartesian_product(E2, E))
        @test X isa EmptySet{N} && dim(X) == 5
    end

    # convex_hull
    E2 = convex_hull(E)
    @test E2 isa EmptySet{N} && dim(E2) == 2
    @test_throws AssertionError convex_hull(E, EmptySet{N}(3))
    E2 = convex_hull(E, E)
    @test E2 isa EmptySet{N} && dim(E2) == 2
    for X in (convex_hull(E, Pnc), convex_hull(Pnc, E))
        @test X isa LazySet{N} && isequivalent(X, Pc)
    end

    # difference
    E2 = EmptySet{N}(3)
    @test_throws AssertionError difference(B, E2)
    @test_throws AssertionError difference(E2, B)
    for X in (difference(E, E), difference(E, B))
        @test X isa EmptySet{N} && dim(X) == 2
    end
    E2 = difference(B, E)
    @test E2 isa BallInf{N} && E2 == B

    # distance
    @test_throws ArgumentError distance(E, E)
    E2 = EmptySet{N}(3)
    @test_throws AssertionError distance(E, E2)
    @test_throws AssertionError distance(E2, E)

    # exact_sum / minkowski_sum
    for f in (exact_sum, minkowski_sum)
        E2 = EmptySet{N}(3)
        @test_throws AssertionError f(E, E2)
        @test_throws AssertionError f(E2, E)
        for X in (f(E, E), f(E, B), f(B, E))
            @test X isa EmptySet{N} && dim(X) == 2
        end
    end

    # intersection
    @test_throws AssertionError convex_hull(E, EmptySet{N}(3))
    for X in (intersection(E, E), intersection(E, B), intersection(B, E))
        @test X isa EmptySet{N} && dim(X) == 2
    end

    # isapprox
    @test E ≈ E
    E2 = EmptySet{N}(3)
    @test_throws AssertionError E ≈ E2
    @test_throws AssertionError E2 ≈ E

    # isdisjoint
    @test_throws AssertionError isdisjoint(E, EmptySet{N}(3))
    @test isdisjoint(E, B) && isdisjoint(B, E) && isdisjoint(E, E)
    for (res, w) in (isdisjoint(E, B, true), isdisjoint(B, E, true), isdisjoint(E, E, true))
        @test res && w isa Vector{N} && w == N[]
    end

    # isequal
    @test E == E
    E2 = EmptySet{N}(3)
    @test !(E == E2) && !(E2 == E)

    # isequivalent
    @test isequivalent(E, E)
    E2 = EmptySet{N}(3)
    @test_throws AssertionError isequivalent(E, E2)
    @test_throws AssertionError isequivalent(E2, E)

    # isstrictsubset
    E2 = EmptySet{N}(3)
    @test_throws AssertionError B ⊂ E2
    @test_throws AssertionError E2 ⊂ B
    for X in (E, B)
        @test !(X ⊂ E)
        res, w = ⊂(X, E, true)
        @test !res && w isa Vector{N} && w == N[]
    end
    @test E ⊂ B
    res, w = ⊂(E, B, true)
    @test res && w isa Vector{N} && length(w) == 2 && w ∉ E && w ∈ B

    # issubset
    E2 = EmptySet{N}(3)
    @test_throws AssertionError B ⊆ E2
    @test_throws AssertionError E2 ⊆ B
    for X in (E, B)
        @test E ⊆ X
        res, w = ⊆(E, X, true)
        @test res && w isa Vector{N} && w == N[]
    end
    @test B ⊈ E
    res, w = ⊆(B, E, true)
    @test !res && w isa Vector{N} && length(w) == 2 && w ∉ E && w ∈ B

    # linear_combination
    @test_throws AssertionError linear_combination(E, EmptySet{N}(3))
    for X in (E, linear_combination(E, Pnc), linear_combination(Pnc, E))
        @test X isa EmptySet{N} && dim(X) == 2
    end

    # minkowski_difference
    E2 = EmptySet{N}(3)
    @test_throws AssertionError minkowski_difference(B, E2)
    @test_throws AssertionError minkowski_difference(E2, B)
    for X in (minkowski_difference(E, E), minkowski_difference(E, B))
        @test X isa EmptySet{N} && dim(X) == 2
    end
    X = minkowski_difference(B, E)
    @test X isa BallInf{N} && X == B

    # chebyshev_center_radius
    @test_throws ArgumentError chebyshev_center_radius(E)
end

for N in [Float64, Float32]
    E = EmptySet{N}(2)

    # is_interior_point
    @test_throws AssertionError is_interior_point(N[0], E)
    @test !is_interior_point(N[0, 0], E)
end
