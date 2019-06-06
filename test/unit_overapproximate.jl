for N in [Float64, Rational{Int}, Float32]
    # overapproximating a set of type T1 with an unsupported type T2 is the
    # identity if T1 = T2
    @test_throws MethodError overapproximate(ZeroSet{N}(2), EmptySet)
    e = EmptySet{N}()
    @test overapproximate(e, EmptySet) == e

    # HPolygon approximation with box directions
    c = N[0, 0]
    b = Ball1(c, N(1))
    p = overapproximate(b, HPolygon)
    for d in [N[1, 0], N[-1, 0]]
        @test σ(d, p)[1] ≈ σ(d, b)[1]
    end
    for d in [N[0, 1], N[0, -1]]
        @test σ(d, p)[2] ≈ σ(d, b)[2]
    end

    # Hyperrectangle approximation
    c = N[0, 0]
    b = Ball1(c, N(1))
    p = overapproximate(b, Hyperrectangle)
    for d in [N[1, 0], N[-1, 0]]
        @test σ(d, p)[1] ≈ σ(d, b)[1]
    end
    for d in [N[0, 1], N[0, -1]]
        @test σ(d, p)[2] ≈ σ(d, b)[2]
    end
    @test p.center ≈ c
    @test p.radius ≈ N[1, 1]

    # Interval approximation
    b = Ball1(N[0], N(1))
    p = overapproximate(b, LazySets.Interval)
    for d in [N[1], N[-1]]
        @test σ(d, p)[1] ≈ σ(d, b)[1]
    end

    # approximation with an axis-aligned hyperrectangle
    Z = Zonotope(N[-1.0, -1.0], N[-1/2 0; -1.2 -1])
    Zoa = overapproximate(Z, Hyperrectangle) # faster o.a.
    Zba = box_approximation(Z) # default o.a. implementation that uses supp function
    @test Zoa.center ≈ Zba.center && Zoa.radius ≈ Zba.radius

    # same but using cartesian product
    h1 = Hyperrectangle(N[1/2],  N[1/2])
    h2 = Hyperrectangle(N[2.5, 4.5],  N[1/2, 1/2])
    H = overapproximate(h1 × h2, Hyperrectangle) # defaults to convert method
    @test low(H) == N[0, 2, 4] && high(H) == N[1, 3, 5]

    # overapproximation of the lazy linear map of a hyperrectangular set
    H = Hyperrectangle(N[0, 0], N[1/2, 1])
    M = Diagonal(N[2, 2])
    OA = overapproximate(M*H, Hyperrectangle)
    @test OA isa Hyperrectangle && OA.center == N[0, 0] && OA.radius == N[1, 2]

    #overapproximation of Minkowski sum of linear maps for each block in the row block
    i1 = Interval(N[0, 1])
    h = Hyperrectangle(low=N[3, 4], high=N[5, 7])
    M = N[1 2 3; 4 5 6; 7 8 9]
    cpa = CartesianProductArray([i1, h])
    lm = M * cpa

    oa = overapproximate(lm, Hyperrectangle)
    oa_box = overapproximate(lm, Approximations.BoxDirections)
    d_oa_d_hp = overapproximate(lm, CartesianProductArray{N, Hyperrectangle{N}})
    d_oa_d_box = overapproximate(lm, CartesianProductArray, Approximations.BoxDirections)
    oa_d_hp = overapproximate(d_oa_d_hp)
    oa_d_box = overapproximate(d_oa_d_box, Approximations.BoxDirections)

    @test oa == oa_d_hp
    @test oa_box == oa_d_box

    for (oax, set_type) in [(d_oa_d_hp, Hyperrectangle), (d_oa_d_box, HPolytope)]
        @test oax isa CartesianProductArray
        arr = oax.array
        @test length(arr) == 2 && dim(arr[1]) == 1 && dim(arr[2]) == 2
        @test all(X -> X isa set_type, arr)
    end

    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    i3 = Interval(N[1, 4])
    cpa = CartesianProductArray([i1, i2, i3])
    M = N[1 2 0; 0 1 0; 0 1 1]
    lm = M * cpa
    d_oa = overapproximate(lm, CartesianProductArray{N, Interval{N}})
    oa = overapproximate(lm)
    @test overapproximate(d_oa) == oa
    @test typeof(d_oa) == CartesianProductArray{N, Interval{N}}

end

# tests that do not work with Rational{Int}
for N in [Float64, Float32]
    # useful for benchmarking overapproximate and LinearMap's support vector
    # (see #290)
    function overapproximate_lmap(n)
        B = BallInf(ones(N, n), N(2))
        π = sparse(N[1, 2], N[1, 2], ones(N, 2), 2, n)
        return Approximations.overapproximate(π*B)
    end
    o = overapproximate_lmap(50)
    @test o.center == N[1, 1] && o.radius == N[2, 2]

    # Approximation of a 2D centered unit ball in norm 1
    # All vertices v should be like this:
    # ‖v‖ >= 1 and ‖v‖ <= 1+ε
    # Where ε is the given error bound
    b = Ball1(N[0, 0], N(1))
    ε = N(0.01)
    p = tovrep(overapproximate(b, ε))
    for v in vertices_list(p)
    @test norm(v) >= N(1)
    @test norm(v) <= N(1 + ε)
    end

    # Check that there are no redundant constraints for a ballinf
    b = BallInf(N[0.5, 0.5], N(0.1))
    lcl = overapproximate(b, N(0.001)).constraints
    @test length(lcl) == 4
    @test lcl[1].a == N[1, 0]
    @test lcl[1].b == N(0.6)
    @test lcl[2].a == N[0, 1]
    @test lcl[2].b == N(0.6)
    @test lcl[3].a == N[-1, 0]
    @test lcl[3].b == N(-0.4)
    @test lcl[4].a == N[0, -1]
    @test lcl[4].b == N(-0.4)

    # Check that there are no redundant constraints for a HPolygon (octagon)
    p = HPolygon{N}()
    addconstraint!(p, LinearConstraint(N[1, 0], N(1)))
    addconstraint!(p, LinearConstraint(N[0, 1], N(1)))
    addconstraint!(p, LinearConstraint(N[-1, 0], N(1)))
    addconstraint!(p, LinearConstraint(N[0, -1], N(1)))
    addconstraint!(p, LinearConstraint(N[sqrt(2)/2, sqrt(2)/2], N(1)))
    addconstraint!(p, LinearConstraint(N[-sqrt(2)/2, sqrt(2)/2], N(1)))
    addconstraint!(p, LinearConstraint(N[sqrt(2)/2, -sqrt(2)/2], N(1)))
    addconstraint!(p, LinearConstraint(N[-sqrt(2)/2, -sqrt(2)/2], N(1)))
    lcl = overapproximate(p, N(.001)).constraints
    @test length(lcl) == 8
    @test lcl[1].a ≈ N[1, 0]
    @test lcl[1].b ≈ N(1)
    @test lcl[2].a ≈ N[sqrt(2)/2, sqrt(2)/2]
    @test lcl[2].b ≈ N(1)
    @test lcl[3].a ≈ N[0, 1]
    @test lcl[3].b ≈ N(1)
    @test lcl[4].a ≈ N[-sqrt(2)/2, sqrt(2)/2]
    @test lcl[4].b ≈ N(1)
    @test lcl[5].a ≈ N[-1, 0]
    @test lcl[5].b ≈ N(1)
    @test lcl[6].a ≈ N[-sqrt(2)/2, -sqrt(2)/2]
    @test lcl[6].b ≈ N(1)
    @test lcl[7].a ≈ N[0, -1]
    @test lcl[7].b ≈ N(1)
    @test lcl[8].a ≈ N[sqrt(2)/2, -sqrt(2)/2]
    @test lcl[8].b ≈ N(1)

    # Zonotope approximation
    Z1 = Zonotope(ones(N, 2), [N[1, 0], N[0, 1], N[1, 1]])
    Z2 = Zonotope(-ones(N, 2), [N[0.5, 1], N[-0.1, 0.9], N[1, 4]])
    Y = ConvexHull(Z1, Z2)
    Y_polygon = overapproximate(Y, N(1e-3)) # overapproximate with a polygon
    Y_zonotope = overapproximate(Y, Zonotope) # overapproximate with a zonotope
    @test Y_polygon ⊆ Y_zonotope
    @test !(Y_zonotope ⊆ Y_polygon)
end

for N in [Float64] # due to sparse vectors: a = sparse(Float32[1 -1; 1 1];); a \ Float32[4, 10]
    #decomposed linear map approximation
    i1 = Interval(N[0, 1])
    i2 = Interval(N[2, 3])
    M = N[1 2; 0 1]
    cpa = CartesianProductArray([i1, i2])
    lm = M * cpa
    d_oa = overapproximate(lm, CartesianProductArray{N, Interval{N}})
    oa = overapproximate(lm, OctDirections)
    @test oa ⊆ d_oa
end
