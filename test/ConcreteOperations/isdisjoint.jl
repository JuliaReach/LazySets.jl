for N in [Float64, Float32, Rational{Int}]
    # using IA types
    X = IA.interval(N(0), N(1)) # IA
    Y = Interval(N(-1), N(2))
    Z = Interval(N(2), N(3))
    res, w = isdisjoint(X, Y, true)
    @test !isdisjoint(X, Y) && !isdisjoint(Y, X) && !res && w ∈ Interval(X) && w ∈ Y
    res, w = isdisjoint(X, Z, true)
    @test isdisjoint(X, Z) && isdisjoint(Z, X) && res && w == N[]

    X = IntervalBox(IA.interval(N(0), N(1)), IA.interval(N(0), N(1)))
    Y = Hyperrectangle(; low=[N(-1), N(-1)], high=[N(2), N(2)])
    @test !isdisjoint(X, Y)
    @test !isdisjoint(Y, X)

    if N == Rational{Int} && test_suite_polyhedra
        for n in [2, 3]
            Z = convert(Zonotope, BallInf(zeros(N, n), N(1)))
            P = convert(HPolyhedron, BallInf(3 * ones(N, n), N(1)))
            res, w = isdisjoint(Z, P, true)
            @test isdisjoint(Z, P) && isdisjoint(P, Z) && res && w == N[]
            P = convert(HPolyhedron, BallInf(ones(N, n), N(1)))
            res, w = isdisjoint(Z, P, true)
            @test !isdisjoint(Z, P) && !isdisjoint(P, Z) && !res && w ∈ Z && w ∈ P
        end
    end
end

for N in [Float64, Float32]
    # case 1: ellipse completely outside
    R1 = Hyperrectangle(N[1.0, -1.0], N[2, 1])
    El1 = Ellipsoid(N[5.0, 5.0], Matrix{N}(I, 2, 2))
    @test isdisjoint(R1, El1)

    # case 2: ellipse completely inside
    R2 = Hyperrectangle(N[0.0, 0.0], N[4.0, 1.0])
    El2 = Ellipsoid(N[1.0, 0.0], Matrix{N}(I, 2, 2))
    @test !isdisjoint(R2, El2)

    # case 3: ellipse overlapping
    R3 = Hyperrectangle(N[0.0, 0.0], N[1.0, 1.0])
    El3 = Ellipsoid(N[2.0, 0.0], Diagonal(N[2.0, 1.0]))
    @test !isdisjoint(R3, El3)

    # case 4: ellipse tangent in one point
    R4 = Hyperrectangle(N[0.0, 0.0], N[1.0, 1.0])
    El4 = Ellipsoid(N[2.0, 0.0], Diagonal(N[1.0, 1.0]))
    @test !isdisjoint(R4, El4)
end