for N in [Float64, Float32, Rational{Int}]
    # Cartesian Product of a centered 1D BallInf and a centered 2D BallInf
    # Here a 3D BallInf
    b1 = BallInf(N[0], N(1))
    b2 = BallInf(N[0, 0], N(1))
    # Test Construction
    cp = CartesianProduct(b1, b2)
    @test cp.X == b1
    @test cp.Y == b2
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

    @test ∈(N[0, 0, 0, 0], cp)
    @test ∈(N[4, 2, 1, 0], cp)
    @test ∈(N[2, 4, -1, 0], cp)
    @test ∈(N[-1, 1, 0.5, 0.7], cp)
    @test ∈(N[2, 3, -0.8, 0.9], cp)
    @test ∈(N[1, 1, -1, 0], cp)
    @test ∈(N[3, 2, 0, 1], cp)
    @test ∈(N[5 / 4, 7 / 4, 1, 1], cp)
    @test !∈(N[4, 1, 0, 0], cp)
    @test !∈(N[5, 2, 0, 0], cp)
    @test !∈(N[3, 4, 0, 0], cp)
    @test !∈(N[-1, 5, 0, 0], cp)
    @test !∈(N[4, 2, 3, 1], cp)
    @test !∈(N[2, 3, 3, 1], cp)
    @test !∈(N[0, 0, 3, 1], cp)
    @test !∈(N[1, 1, 3, 1], cp)

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
    @test ∈(N[1, 2, 3, 4], cpa)
    @test !∈(N[3, 4, 1, 2], cpa)

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

    #constraints_list
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
