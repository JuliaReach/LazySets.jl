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
    @test isbounded(cp)
    @test !isbounded(Singleton(N[1]) * HalfSpace(N[1], N(1)))

    # isempty
    @test !isempty(cp)

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

    # Test Cartesian Product with EmptySet
    s = Singleton(N[1])
    E = EmptySet{N}()
    cs1 = E * s
    cs2 = s * E
    @test cs1 isa EmptySet
    @test cs2 isa EmptySet

    # Test Cartesian Product of an array
    # 0-elements
    as = LazySet{N}[]
    cs = CartesianProduct(as)
    @test cs isa EmptySet
    # 1-element
    as = [Singleton(N[1])]
    cs = CartesianProduct(as)
    @test cs.element == N[1]
    # 3-elements
    as = [Singleton(N[1]), Singleton(N[2]), Singleton(N[3])]
    cs = CartesianProduct(as)
    @test cs.X.element == N[1]
    @test cs.Y.X.element == N[2]
    @test cs.Y.Y.element == N[3]

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

    #constraints_list
    hlist = constraints_list(CartesianProduct(i1, i2))
    @test ispermutation(hlist,
        [LinearConstraint(sparsevec(N[1], N[1], 2), N(1)),
        LinearConstraint(sparsevec(N[1], N[-1], 2), N(0)),
        LinearConstraint(sparsevec(N[2], N[1], 2),  N(3)),
        LinearConstraint(sparsevec(N[2], N[-1], 2), N(-2))])
    @test all(H -> dim(H) == 2, hlist)

    # ==================================
    # Conversions of Cartesian Products
    # ==================================

    # convert the cartesian product of two hyperrectangles to one hyperrectangle
    h1 = Hyperrectangle(N[1/2],  N[1/2])
    h2 = Hyperrectangle(N[2.5, 4.5],  N[1/2, 1/2])
    H = convert(Hyperrectangle, h1 × h2)
    @test low(H) == N[0, 2, 4] && high(H) == N[1, 3, 5]

    # make one of the sets a BallInf, to test that the conversion does
    # not depend on the Hyperrectangle type
    b1 = BallInf(N[1/2],  N(1/2))
    H = convert(Hyperrectangle, b1 × h2)
    @test low(H) == N[0, 2, 4] && high(H) == N[1, 3, 5]

    # convert a hyperrectangle to a cartesian product of intervals
    H = Hyperrectangle(N[0, 0], N[1, 1])
    Hcp = convert(CartesianProduct{N, Interval{N}, Interval{N}}, H)
    @test Hcp isa CartesianProduct &&
          Hcp.X == Interval(N(-1), N(1)) &&
          Hcp.Y == Interval(N(-1), N(1))

    # checking that the order is correct
    H = Hyperrectangle(N[0, 1], N[1, 3])
    Hcp = convert(CartesianProduct{N, Interval{N}, Interval{N}}, H)
    @test Hcp isa CartesianProduct &&
          Hcp.X == Interval(N(-1), N(1)) &&
          Hcp.Y == Interval(N(-2), N(4))

    # =====================
    # CartesianProductArray
    # =====================

    # relation to base type (internal helper functions)
    @test LazySets.array_constructor(CartesianProduct) == CartesianProductArray
    @test LazySets.is_array_constructor(CartesianProductArray)

    # standard constructor
    v = Vector{LazySet{N}}()
    push!(v, Singleton(N[1, 2]))
    push!(v, Singleton(N[3, 4]))
    cpa = CartesianProductArray(v)

    # constructor with size hint and type
    CartesianProductArray(10, N)

    # array getter
    @test array(cpa) ≡ v

    # boundedness
    @test isbounded(cpa)
    @test !isbounded(CartesianProductArray([Singleton(N[1]), HalfSpace(N[1], N(1))]))

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
    @test ispermutation(vlist, [N[0, 2, 4], N[0, 3, 4], N[1, 2, 4],
        N[1, 3, 4], N[0, 2, 5], N[0, 3, 5], N[1, 2, 5], N[1, 3, 5]])

    # constraints_list
    hlist = constraints_list(CartesianProductArray([i1, i2, i3]))
    @test ispermutation(hlist,
        [LinearConstraint(sparsevec(N[1], N[1], 3), N(1)),
        LinearConstraint(sparsevec(N[1], N[-1], 3), N(0)),
        LinearConstraint(sparsevec(N[2], N[1], 3), N(3)),
        LinearConstraint(sparsevec(N[2], N[-1], 3), N(-2)),
        LinearConstraint(sparsevec(N[3], N[1], 3), N(5)),
        LinearConstraint(sparsevec(N[3], N[-1], 3), N(-4)),
        ])
    @test all(H -> dim(H) == 3, hlist)
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
    h1 = Hyperrectangle(N[1/2],  N[1/2])
    h2 = Hyperrectangle(N[2.5, 4.5],  N[1/2, 1/2])
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

    # convert a hyperrectangle to the cartesian product array of intervals
    # convert a hyperrectangle to a cartesian product of intervals
    H = Hyperrectangle(N[0, 0, 0], N[1, 1, 1])
    Hcp = convert(CartesianProductArray{N, Interval{N}}, H)
    @test Hcp isa CartesianProductArray &&
          all([Hcpi == Interval(N(-1), N(1)) for Hcpi in array(Hcp)])

    # check that the order is correct
    H = Hyperrectangle(N[0, 1, 2], N[1, 2, 3])
    Hcp = convert(CartesianProductArray{N, Interval{N}}, H)
    @test Hcp isa CartesianProductArray &&
          all([array(Hcp)[i] == Interval(N(-1), N(2*i-1)) for i in 1:3])

    # ================
    # common functions
    # ================

    # absorbing element
    e = EmptySet{N}()
    b = BallInf(N[0, 0], N(2))
    @test absorbing(CartesianProduct) == absorbing(CartesianProductArray) ==
          EmptySet
    @test b × e == e × b == cpa × e == e × cpa == e × e == e
end

for N in [Float64, Float32]
    # is_intersection_empty
    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    h1 = Hyperrectangle(low=N[3, 4], high=N[5, 7])
    h2 = Hyperrectangle(low=N[5, 5], high=N[6, 8])
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
end
