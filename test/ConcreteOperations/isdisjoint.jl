for N in @tN([Float64, Float32, Rational{Int}])
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
