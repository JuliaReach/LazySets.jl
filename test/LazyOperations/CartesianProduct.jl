using LazySets.Approximations: overapproximate

for N in [Float64, Float32, Rational{Int}]
    # Cartesian Product of a centered 1D BallInf and a centered 2D BallInf
    # Here a 3D BallInf
    b1 = BallInf(N[0], N(1))
    b2 = BallInf(N[0, 0], N(1))
    # Test Construction
    cp = CartesianProduct(b1, b2)
    @test cp.X == b1
    @test cp.Y == b2

    # convenience constructors
    @test b1 * b2 == b1 × b2 == cp
    cpa = CartesianProductArray([b1, b2, b1])
    @test b1 * b2 * b1 == *(b1, b2, b1) == ×(b1, b2, b1) == *([b1, b2, b1]) == ×([b1, b2, b1]) ==
          cpa
    @test *(b1) == ×(b1) == b1

    # array interface
    @test array(cp) == [b1, b2] && array(cpa) == [b1, b2, b1]
    @test cp[1] == cpa[1] == b1
    @test length(cp) == 2 && length(cpa) == 3
    v = Vector{LazySet{N}}()
    @test array(CartesianProductArray(v)) ≡ v

    # swap
    cp2 = swap(cp)
    @test cp.X == cp2.Y && cp.Y == cp2.X

    # Test Dimension
    @test dim(cp) == 3
    # Test Support Vector
    d = N[1, 1, 1]
    @test σ(d, cp) == N[1, 1, 1]
    d = N[1, 1, -1]
    @test σ(d, cp) == N[1, 1, -1]
    d = N[1, -1, 1]
    @test σ(d, cp) == N[1, -1, 1]
    d = N[1, -1, -1]
    @test σ(d, cp) == N[1, -1, -1]
    d = N[-1, 1, 1]
    @test σ(d, cp) == N[-1, 1, 1]
    d = N[-1, 1, -1]
    @test σ(d, cp) == N[-1, 1, -1]
    d = N[-1, -1, 1]
    @test σ(d, cp) == N[-1, -1, 1]
    d = N[-1, -1, -1]
    @test σ(d, cp) == N[-1, -1, -1]

    # boundedness
    @test isbounded(cp) && isboundedtype(typeof(cp))
    cp2 = Singleton(N[1]) * HalfSpace(N[1], N(1))
    @test !isbounded(cp2) && !isboundedtype(typeof(cp2))

    # is_polyhedral
    @test is_polyhedral(cp)
    if N isa AbstractFloat
        cp2 = CartesianProduct(b1, Ball2(N[0, 0], N(1)))
        @test !is_polyhedral(cp2)
    end

    # isempty
    @test !isempty(cp)

    # center
    @test center(cp) == N[0, 0, 0]

    # Cartesian Product of a not-centered 1D BallInf and a not-centered 2D BallInf
    # Here a Hyperrectangle where c = [1, -3, 4] and r = [3, 2, 2]
    b1 = BallInf(N[1], N(3))
    b2 = BallInf(N[-3, 4], N(2))
    # Test Construction
    c = CartesianProduct(b1, b2)
    # Test Dimension
    @test dim(c) == 3
    # Test Support Vector
    d = N[1, 1, 1]
    @test σ(d, c) == N[4, -1, 6]
    d = N[1, 1, -1]
    @test σ(d, c) == N[4, -1, 2]
    d = N[1, -1, 1]
    @test σ(d, c) == N[4, -5, 6]
    d = N[1, -1, -1]
    @test σ(d, c) == N[4, -5, 2]
    d = N[-1, 1, 1]
    @test σ(d, c) == N[-2, -1, 6]
    d = N[-1, 1, -1]
    @test σ(d, c) == N[-2, -1, 2]
    d = N[-1, -1, 1]
    @test σ(d, c) == N[-2, -5, 6]
    d = N[-1, -1, -1]
    @test σ(d, c) == N[-2, -5, 2]

    # Test containment with respect to CartesianProduct
    p1 = HPolygon{N}()
    addconstraint!(p1, LinearConstraint(N[2, 2], N(12)))
    addconstraint!(p1, LinearConstraint(N[-3, 3], N(6)))
    addconstraint!(p1, LinearConstraint(N[-1, -1], N(0)))
    addconstraint!(p1, LinearConstraint(N[2, -4], N(0)))
    p2 = HPolygon{N}()
    addconstraint!(p2, LinearConstraint(N[1, 0], N(1)))
    addconstraint!(p2, LinearConstraint(N[-1, 0], N(1)))
    addconstraint!(p2, LinearConstraint(N[0, 1], N(1)))
    addconstraint!(p2, LinearConstraint(N[0, -1], N(1)))
    cp = CartesianProduct(p1, p2)

    @test N[0, 0, 0, 0] ∈ cp
    @test N[4, 2, 1, 0] ∈ cp
    @test N[2, 4, -1, 0] ∈ cp
    @test N[-1, 1, 0.5, 0.7] ∈ cp
    @test N[2, 3, -0.8, 0.9] ∈ cp
    @test N[1, 1, -1, 0] ∈ cp
    @test N[3, 2, 0, 1] ∈ cp
    @test N[5 / 4, 7 / 4, 1, 1] ∈ cp
    @test N[4, 1, 0, 0] ∉ cp
    @test N[5, 2, 0, 0] ∉ cp
    @test N[3, 4, 0, 0] ∉ cp
    @test N[-1, 5, 0, 0] ∉ cp
    @test N[4, 2, 3, 1] ∉ cp
    @test N[2, 3, 3, 1] ∉ cp
    @test N[0, 0, 3, 1] ∉ cp
    @test N[1, 1, 3, 1] ∉ cp

    # vertices_list
    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    vlist = vertices_list(CartesianProduct(i1, i2))
    @test ispermutation(vlist, [N[0, 2], N[0, 3], N[1, 2], N[1, 3]])

    # constraints_list
    hlist = constraints_list(CartesianProduct(i1, i2))
    @test ispermutation(hlist,
                        [LinearConstraint(sparsevec(N[1], N[1], 2), N(1)),
                         LinearConstraint(sparsevec(N[1], N[-1], 2), N(0)),
                         LinearConstraint(sparsevec(N[2], N[1], 2), N(3)),
                         LinearConstraint(sparsevec(N[2], N[-1], 2), N(-2))])
    @test all(H -> dim(H) == 2, hlist)
    # empty constraints list is handled correctly
    H = HalfSpace(N[1], N(0))
    hlist = constraints_list(CartesianProduct(H, Universe{N}(2)))
    @test hlist == [HalfSpace(N[1, 0, 0], N(0))]

    # linear_map
    cp = CartesianProduct(Interval(N(0), N(1)), Interval(N(2), N(3)))
    M = N[2 1; 0 2]  # invertible matrix
    lm = linear_map(M, cp)
    @test lm isa HPolytope{N} && length(constraints_list(lm)) ==
                                 length(constraints_list(cp)) == 4
    M = N[2 1; 0 0]  # singular matrix
    lm = linear_map(M, cp)
    @test lm isa (N == Float64 ? HPolytope{N} : HPolytope)
    @test box_approximation(lm) == Hyperrectangle(N[7 // 2, 0], N[3 // 2, 0])

    cp = CartesianProduct(VPolytope([N[1]]), VPolytope([N[2]]))
    if test_suite_polyhedra
        # concretize
        @test concretize(cp) == VPolytope([N[1, 2]])
    else
        @test concretize(cp) === cp
    end

    # inclusion
    # same block structure
    cp1 = CartesianProduct(Ball1(N[2], N(1)), Ball1(N[0, 0], N(2)))
    cp2 = CartesianProduct(Interval(N(1), N(3)), BallInf(N[0, 0], N(3)))
    res, w = ⊆(cp1, cp2, true)
    @test cp1 ⊆ cp2 && res && w == N[]
    res, w = ⊆(cp2, cp1, true)
    @test cp2 ⊈ cp1 && !res && w ∈ cp2 && w ∉ cp1
    # different block structure
    cp1 = CartesianProduct(Interval(N(1), N(2)), BallInf(N[1, 1], N(1)))
    cp2 = CartesianProduct(BallInf(N[1, 1], N(2)), Interval(N(0), N(2)))
    res, w = ⊆(cp1, cp2, true)
    @test cp1 ⊆ cp2 && res && w == N[]
    res, w = ⊆(cp2, cp1, true)
    @test cp2 ⊈ cp1 && !res && w ∈ cp2 && w ∉ cp1

    # projection
    P = Singleton(N[11, 12])
    Q = Singleton(N[13, 14, 15])
    cp = P × Q
    @test project(cp, [2]) == Singleton(N[12])
    @test project(cp, [3, 5]) == Singleton(N[13, 15])
    @test project(cp, [1, 4, 5]) == Singleton(N[11]) × Singleton(N[14, 15])

    # ==================================
    # Conversions of Cartesian Products
    # ==================================

    # convert the cartesian product of two hyperrectangles to one hyperrectangle
    h1 = Hyperrectangle(N[1 / 2], N[1 / 2])
    h2 = Hyperrectangle(N[2.5, 4.5], N[1 / 2, 1 / 2])
    H = convert(Hyperrectangle, h1 × h2)
    @test low(H) == N[0, 2, 4] && high(H) == N[1, 3, 5]

    # make one of the sets a BallInf, to test that the conversion does
    # not depend on the Hyperrectangle type
    b1 = BallInf(N[1 / 2], N(1 / 2))
    H = convert(Hyperrectangle, b1 × h2)
    @test low(H) == N[0, 2, 4] && high(H) == N[1, 3, 5]

    # convert a hyperrectangle to a cartesian product of intervals
    H = Hyperrectangle(N[0, 0], N[1, 1])
    Hcp = convert(CartesianProduct{N,Interval{N},Interval{N}}, H)
    @test Hcp isa CartesianProduct &&
          Hcp.X == Interval(N(-1), N(1)) &&
          Hcp.Y == Interval(N(-1), N(1))

    # checking that the order is correct
    H = Hyperrectangle(N[0, 1], N[1, 3])
    Hcp = convert(CartesianProduct{N,Interval{N},Interval{N}}, H)
    @test Hcp isa CartesianProduct &&
          Hcp.X == Interval(N(-1), N(1)) &&
          Hcp.Y == Interval(N(-2), N(4))

    cp = Zonotope(N[0, 0], N[1 0; 0 1]) × Zonotope(N[1, 1], N[0 1; 1 0])
    @test convert(Zonotope, cp) == Zonotope(N[0, 0, 1, 1], N[1 0 0 0; 0 1 0 0;
                                                             0 0 0 1; 0 0 1 0])

    # cartesian product of singletons
    S1 = Singleton(N[1, 2, 3])
    S2 = Singleton(N[4, 5, 6])
    @test convert(Singleton, S1 × S2) == Singleton(N[1, 2, 3, 4, 5, 6])

    # =====================
    # CartesianProductArray
    # =====================

    # relation to base type (internal helper functions)
    @test LazySets.array_constructor(CartesianProduct) == CartesianProductArray
    @test LazySets.binary_constructor(CartesianProductArray) == CartesianProduct
    @test !LazySets.is_array_constructor(CartesianProduct)
    @test LazySets.is_array_constructor(CartesianProductArray)

    # standard constructor
    v = Vector{Singleton{N}}()
    S1 = Singleton(N[1, 2])
    S2 = Singleton(N[3, 4])
    push!(v, S1)
    push!(v, S2)
    cpa = CartesianProductArray(v)

    # constructor with size hint and type
    CartesianProductArray(10, N)

    # getindex & length
    @test cpa[1] == S1 && cpa[2] == S2
    @test length(cpa) == 2

    # support function & support vector
    svec = N[1, 2, 3, 4]
    for dd in [N[1, 1, 0, 0], N[0, 0, 1, 1], N[0, 1, 0, 1], N[0, 0, 1, 0], N[0, 0, 0, 0]]
        ds = sparse(dd)
        @test σ(dd, cpa) == σ(ds, cpa) == svec
        @test ρ(dd, cpa) == ρ(ds, cpa) == dot(svec, dd)
    end

    # center
    @test center(cpa) == N[1, 2, 3, 4]

    # boundedness
    @test isbounded(cpa) && isboundedtype(typeof(cpa))
    cpa2 = CartesianProductArray([Singleton(N[1]), HalfSpace(N[1], N(1))])
    @test !isbounded(cpa2) && !isboundedtype(typeof(cpa2))

    # membership
    @test N[1, 2, 3, 4] ∈ cpa
    @test N[3, 4, 1, 2] ∉ cpa

    # in-place modification
    b = BallInf(N[0, 0], N(2))
    cpa = CartesianProductArray(LazySet{N}[])
    @test CartesianProduct!(b, b) isa CartesianProduct &&
          length(array(cpa)) == 0
    res = CartesianProduct!(b, cpa)
    @test res isa CartesianProductArray && length(array(cpa)) == 1
    res = CartesianProduct!(cpa, b)
    @test res isa CartesianProductArray && length(array(cpa)) == 2
    res = CartesianProduct!(cpa, cpa)
    @test res isa CartesianProductArray && length(array(cpa)) == 4

    # vertices_list
    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    i3 = Interval(N[4, 5])
    vlist = vertices_list(CartesianProductArray([i1, i2, i3]))
    @test ispermutation(vlist,
                        [N[0, 2, 4], N[0, 3, 4], N[1, 2, 4],
                         N[1, 3, 4], N[0, 2, 5], N[0, 3, 5], N[1, 2, 5], N[1, 3, 5]])

    # constraints_list
    hlist = constraints_list(CartesianProductArray([i1, i2, i3]))
    @test ispermutation(hlist,
                        [LinearConstraint(sparsevec(N[1], N[1], 3), N(1)),
                         LinearConstraint(sparsevec(N[1], N[-1], 3), N(0)),
                         LinearConstraint(sparsevec(N[2], N[1], 3), N(3)),
                         LinearConstraint(sparsevec(N[2], N[-1], 3), N(-2)),
                         LinearConstraint(sparsevec(N[3], N[1], 3), N(5)),
                         LinearConstraint(sparsevec(N[3], N[-1], 3), N(-4))])
    @test all(H -> dim(H) == 3, hlist)
    # empty constraints list is handled correctly
    H = HalfSpace(N[1], N(0))
    hlist = constraints_list(CartesianProductArray([H, Universe{N}(2)]))
    @test hlist == [HalfSpace(N[1, 0, 0], N(0))]

    # is_intersection_empty
    cpa = CartesianProductArray([BallInf(ones(N, 5), N(1)),
                                 BallInf(2 * ones(N, 5), N(0.5))])
    b = BallInf(zeros(N, 10), N(1))
    for (x, y) in [(cpa, b), (b, cpa)]
        res, w = isdisjoint(x, y, true)
        @test isdisjoint(x, y) && res && w == N[]
    end
    b = BallInf(ones(N, 10), N(1))
    for (x, y) in [(cpa, b), (b, cpa)]
        res, w = isdisjoint(x, y, true)
        @test !isdisjoint(x, y) && !res && w ∈ x && w ∈ y
    end

    # concrete intersection where all dimensions are unconstrained
    @test intersection(cpa, HPolyhedron{N}()) == cpa

    # linear_map
    cpa = CartesianProductArray([Interval(N(0), N(1)), Interval(N(2), N(3))])
    M = N[2 1; 0 2]  # invertible matrix
    lm = linear_map(M, cpa)
    @test lm isa HPolytope{N} && length(constraints_list(lm)) ==
                                 length(constraints_list(cpa)) == 4
    M = N[2 1; 0 0]  # singular matrix
    lm = linear_map(M, cpa)
    @test lm isa (N == Float64 ? HPolytope{N} : HPolytope)
    @test box_approximation(lm) == Hyperrectangle(N[7 // 2, 0], N[3 // 2, 0])

    cpa = CartesianProductArray([VPolytope([N[1]]), VPolytope([N[2]])])
    if test_suite_polyhedra
        # concretize
        @test concretize(cpa) == VPolytope([N[1, 2]])
    else
        @test concretize(cpa) === cpa
    end

    # projection
    P = Singleton(N[11, 12])
    Q = Singleton(N[13, 14, 15])
    R = Singleton(N[16, 17])
    cpa = CartesianProductArray([P, Q, R])
    @test project(cpa, [2]) == Singleton(N[12])
    @test project(cpa, [4]) == Singleton(N[14])
    @test project(cpa, [6]) == Singleton(N[16])
    @test project(cpa, [3, 5]) == Singleton(N[13, 15])
    @test project(cpa, [2, 7]) == Singleton(N[12]) × Singleton(N[17])
    @test project(cpa, [1, 4, 5]) == Singleton(N[11]) × Singleton(N[14, 15])
    @test project(cpa, [1, 4, 7]) ==
          CartesianProductArray([Singleton(N[11]), Singleton(N[14]), Singleton(N[17])])

    # ========================================
    # Conversions of Cartesian Product Arrays
    # ========================================

    # conversion of the cartesian product array of intervals to a hyperrectangle
    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    i3 = Interval(N[4, 5])
    H = convert(Hyperrectangle, CartesianProductArray([i1, i2, i3]))
    @test low(H) == N[0, 2, 4] && high(H) == N[1, 3, 5]

    M = hcat(N(1))
    h1 = Hyperrectangle(N[1 / 2], N[1 / 2])
    h2 = Hyperrectangle(N[2.5, 4.5], N[1 / 2, 1 / 2])
    H = convert(Hyperrectangle, CartesianProductArray([h1, h2]))
    @test low(H) == N[0, 2, 4] && high(H) == N[1, 3, 5]

    # same block structure
    aX = [Ball1(N[2], N(1)), Ball1(N[0, 0], N(2))]
    aY = [Interval(N(1), N(3)), BallInf(N[0, 0], N(3))]
    @test LazySets.same_block_structure(aX, aY)

    # inclusion
    cpa1 = CartesianProductArray(aX)
    cpa2 = CartesianProductArray(aY)
    res, w = ⊆(cpa1, cpa2, true)
    @test cpa1 ⊆ cpa2 && res && w == N[]
    res, w = ⊆(cpa2, cpa1, true)
    @test cpa2 ⊈ cpa1 && !res && w ∈ cpa2 && w ∉ cpa1
    # different block structure
    cpa1 = CartesianProductArray([Interval(N(1), N(2)), BallInf(N[1, 1], N(1))])
    cpa2 = CartesianProductArray([BallInf(N[1, 1], N(2)), Interval(N(0), N(2))])
    res, w = ⊆(cpa1, cpa2, true)
    @test cpa1 ⊆ cpa2 && res && w == N[]
    res, w = ⊆(cpa2, cpa1, true)
    @test cpa2 ⊈ cpa1 && !res && w ∈ cpa2 && w ∉ cpa1

    # convert a hyperrectangle to the cartesian product array of intervals
    # convert a hyperrectangle to a cartesian product of intervals
    H = Hyperrectangle(N[0, 0, 0], N[1, 1, 1])
    Hcp = convert(CartesianProductArray{N,Interval{N}}, H)
    @test Hcp isa CartesianProductArray &&
          all([Hcpi == Interval(N(-1), N(1)) for Hcpi in array(Hcp)])

    # check that the order is correct
    H = Hyperrectangle(N[0, 1, 2], N[1, 2, 3])
    Hcp = convert(CartesianProductArray{N,Interval{N}}, H)
    @test Hcp isa CartesianProductArray &&
          all([array(Hcp)[i] == Interval(N(-1), N(2 * i - 1)) for i in 1:3])

    cpa = CartesianProductArray([Zonotope(N[0, 0], N[1 0; 0 1]),
                                 Zonotope(N[1, 1], N[0 1; 1 0])])
    @test convert(Zonotope, cpa) == Zonotope(N[0, 0, 1, 1], N[1 0 0 0; 0 1 0 0;
                                                              0 0 0 1; 0 0 1 0])

    # ================
    # common functions
    # ================

    # absorbing element
    e = EmptySet{N}(2)
    b = BallInf(N[0, 0], N(2))
    @test b × e == e × b == e × e == EmptySet{N}(4)
    @test cpa × e == e × cpa == EmptySet{N}(6)

    # volume
    b2 = BallInf(zeros(N, 3), N(3))
    cp = CartesianProduct(b, b2)
    cpa = CartesianProductArray([b, b2])
    @test volume(cp) == volume(cpa) == N(3456)
end

for N in [Float64, Float32]
    # is_intersection_empty
    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    h1 = Hyperrectangle(; low=N[3, 4], high=N[5, 7])
    h2 = Hyperrectangle(; low=N[5, 5], high=N[6, 8])
    cpa1 = CartesianProductArray([i1, i2, h1])
    cpa2 = CartesianProductArray([i1, i2, h2])
    G = HalfSpace(N[1, 0, 0, 0], N(1))
    G_empty = HalfSpace(N[1, 0, 0, 0], N(-1))
    cpa1_box = overapproximate(cpa1)
    cpa2_box = overapproximate(cpa2)

    @test !is_intersection_empty(cpa1, cpa2) &&
          !is_intersection_empty(cpa1_box, cpa2_box)
    @test is_intersection_empty(cpa1, G_empty) &&
          is_intersection_empty(cpa1_box, G_empty)
    @test !is_intersection_empty(cpa1, G) &&
          !is_intersection_empty(Approximations.overapproximate(cpa1), G)
    @test !is_intersection_empty(cpa1, Universe{N}(4)) &&
          !is_intersection_empty(Universe{N}(4), cpa2)
    @test_throws AssertionError is_intersection_empty(cpa1, Universe{N}(3))
    @test_throws AssertionError is_intersection_empty(Universe{N}(5), cpa2)

    # projection
    cp = Interval(N(0), N(2)) × Hyperrectangle(N[2, 3], N[1, 1])
    @test project(cp, 1:2) == Hyperrectangle(N[1, 2], N[1, 1])
    @test project(cp, 2:3) == Hyperrectangle(N[2, 3], N[1, 1])

    cp = Interval(N(0), N(2)) × Zonotope(N[2, 3], N[1 0; 0 1])
    @test isequivalent(project(cp, 1:2), Zonotope(N[1, 2], N[1 0; 0 1]))
    @test isequivalent(project(cp, 2:3), Zonotope(N[2, 3], N[1 0; 0 1]))
end

for N in [Float64]
    # concrete intersection
    # constrained dimensions cover all dimensions
    cpa = CartesianProductArray([Interval(N[0, 1]), Interval(N[1, 2]),
                                 Interval(N[2, 3])])
    P = HalfSpace(N[1, 0, 1], N(2))
    cap = intersection(cpa, P)
    @test cap isa CartesianProductArray && length(array(cap)) == 1
    Q = array(cap)[1]
    @test ispermutation(constraints_list(Q),
                        [HalfSpace(N[-1, 0, 0], N(0)),
                         HalfSpace(N[0, 1, 0], N(2)), HalfSpace(N[0, -1, 0], N(-1)),
                         HalfSpace(N[0, 0, -1], N(-2)), HalfSpace(N[1, 0, 1], N(2))])
    # some dimensions are unconstrained
    cpa = CartesianProductArray([Interval(N[0, 1]),
                                 Hyperrectangle(N[2, 3], N[1, 1]),
                                 Interval(N[3, 4])])
    P = HalfSpace(N[0, 1, 1, 0], N(3))
    cap = intersection(cpa, P)
    @test cap isa CartesianProductArray && length(array(cap)) == 3
    for i in [1, 3]
        @test array(cap)[1] === array(cpa)[1]
    end
    Q = array(cap)[2]
    @test ispermutation(constraints_list(Q),
                        [HalfSpace(N[-1, 0], N(-1)),
                         HalfSpace(N[0, -1], N(-2)), HalfSpace(N[1, 1], N(3))])

    if test_suite_polyhedra
        # projection to mixed dimensions
        P = Singleton(N[11, 12])
        Q = Singleton(N[13, 14, 15])
        cp = P × Q
        @test isequivalent(project(cp, [1, 4]), Singleton(N[11, 14]))
    end
end
