using Plots, Optim

for N in [Float64, Rational{Int}, Float32]
    p0 = zero(N)
    p1 = one(N)
    v0 = zeros(N, 2)
    v1 = ones(N, 2)

    # ---------------
    # set definitions
    # ---------------

    # bounded basic set types
    b1 = Ball1(v0, p1)
    bi = BallInf(v0, p1)
    hr = Hyperrectangle(v0, v1)
    itv = Interval(p0, p1)
    ls = LineSegment(v0, v1)
    st = Singleton(v1)
    zt = Zonotope(v0, Diagonal(N[1, 1]))
    zs = ZeroSet{N}(2)

    # polygon/polytope types
    constraints = [LinearConstraint([p1, p0], p1),
                   LinearConstraint([p0, p1], p1),
                   LinearConstraint([-p1, p0], p0),
                   LinearConstraint([p0, -p1], p0)]
    hpg = HPolygon(constraints)
    hpgo = HPolygonOpt(constraints)
    hpt = HPolytope(constraints)
    vertices = vertices_list(bi)
    vpg = VPolygon(vertices)
    vpt = VPolytope(vertices)

    # empty set
    es = EmptySet{N}()

    # infinite sets
    hs = HalfSpace(v1, p1)
    hp = Hyperplane(v1, p1)
    l = Line(v1, p1)

    # unary set operations
    sih = SymmetricIntervalHull(b1)
    lm = LinearMap(N[2 1; 1 2], bi)

    # binary set operations
    ch = ConvexHull(b1, bi)
    cha = ConvexHullArray([b1, bi])

    ms = MinkowskiSum(b1, bi)
    msa = MinkowskiSumArray([b1, bi])
    cms = CacheMinkowskiSum([b1, bi])
    cp = CartesianProduct(itv, itv)
    cpa = CartesianProductArray([itv, itv])

    # -----
    # plots
    # -----

    # direct
    plot(b1)
    plot(bi)
    plot(hr)
    plot(itv)
    plot(ls)
    plot(st)
    plot(zt)
    plot(zs)
    plot(hpg)
    plot(hpgo)
    plot(hpt)
    plot(vpg)
    plot(vpt)
    @test_throws AssertionError plot(es) # TODO see #577
    @test_throws Exception plot(hs) # TODO see #576
    @test_throws Exception plot(hp) # TODO see #576
    @test_throws ErrorException plot(l) # TODO see #576
    plot(ch)
    plot(cha)
    plot(sih)
    plot(lm)
    plot(ms)
    plot(msa)
    plot(cms)
    plot(cp)
    plot(cpa)

    # ε-close
    ε = N(1e-2)
    if N == Float64 # TODO see #578
        plot(b1, ε)
        plot(bi, ε)
        plot(hr, ε)
        plot(ls, ε)
        @test_throws AssertionError plot(itv, ε) # TODO see #575
        plot(st, ε)
        plot(zt, ε)
        plot(zs, ε)
        plot(hpg, ε)
        plot(hpgo, ε)
        plot(hpt, ε)
        plot(vpg, ε)
        plot(vpt, ε)
        @test_throws AssertionError plot(es, ε) # TODO see #577
        @test_throws ErrorException plot(hs, ε) # TODO see #576/#578
        @test_throws ErrorException plot(hp, ε) # TODO see #576/#578
        @test_throws ErrorException plot(l, ε) # TODO see #576/#578
        plot(ch, ε)
        plot(cha, ε)
        plot(sih, ε)
        plot(lm, ε)
        plot(ms, ε)
        plot(msa, ε)
        plot(cms, ε)
        plot(cp, ε)
        plot(cpa, ε)
    else
        @test_throws ErrorException plot(b1, ε) # TODO see #578
    end
end

# set types that do not work with Rational{Int}
for N in [Float64, Float32]
    v0 = zeros(N, 2)
    p1 = one(N)

    # ---------------
    # set definitions
    # ---------------

    # bounded basic set types
    b2 = Ball2(v0, p1)
    bp = Ballp(N(1.5), v0, p1)
    el = Ellipsoid(v0, Diagonal(N[1, 1]))

    # unary set operations
    spI = SparseMatrixCSC{N}(2I, 2, 2)
    sme = SparseMatrixExp(spI)
    em = ExponentialMap(sme, b2)
    psme = ProjectionSparseMatrixExp(spI, sme, spI)
    epm = ExponentialProjectionMap(psme, b2)

    # -----
    # plots
    # -----

    # direct
    plot(b2)
    plot(bp)
    plot(el)
    plot(em)
    plot(epm)

    # ε-close
    ε = N(1e-2)
    if N == Float64 # TODO see #578
        plot(b2, ε)
        plot(bp, ε)
        plot(el, ε)
        plot(em, ε)
        plot(epm, ε)
    else
        @test_throws ErrorException plot(b2, ε) # TODO see #578
    end
end

# set types that only work with Float64
for N in [Float64]
    v0 = zeros(N, 2)
    p1 = one(N)

    b1 = Ball1(v0, p1)
    bi = BallInf(v0, p1)

    # binary set operations
    its = Intersection(b1, bi)
    itsa = IntersectionArray([b1, bi])

    # -----
    # plots
    # -----

    plot(its)
    @test_throws ErrorException plot(itsa) # TODO not implemented yet

    # ε-close
    ε = N(1e-2)
    plot(its)
    @test_throws ErrorException plot(itsa, ε) # TODO not implemented yet
end
