using SparseArrays

for N in [Float64, Float32, Rational{Int}]
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

    # -------------------------------------------
    # plot polytopes
    # -------------------------------------------
    plot(b1)
    plot(bi)
    plot(hr)
    plot(itv)
    plot(ls)
    plot(st)
    plot(zt)
    plot(zs)
    
    if N == Rational{Int}
        # rationals do not support epsilon-close approximation
        continue
    end

    # bounded basic set types which are not polytopes
    b2 = Ball2(v0, p1)
    bp = Ballp(N(1.5), v0, p1)
    el = Ellipsoid(v0, Diagonal(N[1, 1]))

    # unary set operations
    spI = SparseMatrixCSC{N}(2I, 2, 2)
    sme = SparseMatrixExp(spI)
    em = ExponentialMap(sme, b2)
    psme = ProjectionSparseMatrixExp(spI, sme, spI)
    epm = ExponentialProjectionMap(psme, b2)

    # binary set operations
    its = Intersection(b1, bi)
    itsa = IntersectionArray([b1, bi])

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
    uni = Universe{N}(2)

    # unary set operations
    sih = SymmetricIntervalHull(b1)
    lm = LinearMap(N[2 1; 1 2], bi)
    rm = ResetMap(bi, Dict(1 => N(1)))

    # binary set operations
    ch = ConvexHull(b1, bi)
    cha = ConvexHullArray([b1, bi])

    ms = MinkowskiSum(b1, bi)
    msa = MinkowskiSumArray([b1, bi])
    cms = CacheMinkowskiSum([b1, bi])
    cp = CartesianProduct(itv, itv)
    cpa = CartesianProductArray([itv, itv])

    # ------------------------------------------------------------------
    # plot using epsilon-close approximation (default threshold ε value)
    # ------------------------------------------------------------------
    if N == Float64 # Float32 requires promotion see #1304 
        plot(its)
    end
    plot(hpg)
    plot(hpgo)
    if test_suite_polyhedra
        plot(hpt)
    end
    plot(vpg)
    plot(vpt)
    plot(es)
    @test_throws ErrorException plot(hs) # TODO see #576
    @test_throws ErrorException plot(hp) # TODO see #576
    @test_throws ErrorException plot(l) # TODO see #576
    @test_throws ErrorException plot(uni) # TODO see #576
    plot(ch)
    plot(cha)
    plot(sih)
    plot(lm)
    plot(rm)
    plot(ms)
    plot(msa)
    plot(cms)
    plot(cp)
    plot(cpa)

    # -----
    # cases that are not implemented
    # -----
    @test_throws ErrorException plot(itsa) # TODO not implemented yet
end
