using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # Zonotope and Hyperplane
    Z = Zonotope(N[1, 1], N[1 1; -1 1])
    # disjoint
    H = Hyperplane(N[1, 1], N(3))
    @test !isdisjoint(Z, H) && !isdisjoint(H, Z)
    for (res, w) in (isdisjoint(Z, H, true), isdisjoint(H, Z, true))
        @test !res && w ∈ Z && w ∈ H
    end
    # overlapping
    H = Hyperplane(N[1, 1], N(5))
    @test isdisjoint(Z, H) && isdisjoint(H, Z)
    for (res, w) in (isdisjoint(Z, H, true), isdisjoint(H, Z, true))
        @test res && w isa Vector{N} && isempty(w)
    end
    # zonotope without generators (#2204)
    Z = Zonotope(N[0, 0], Matrix{N}(undef, 2, 0))
    @test isdisjoint(Z, H)

    # Zonotope and HPolyhedron
    @static if isdefined(@__MODULE__, :Polyhedra)
        if N == Rational{Int}
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

    # AbstractSingleton and LazySet
    S = Singleton(N[2, -1])
    H = Hyperrectangle(N[2, 0], N[1, 1])
    L = LinearMap(N[1 0; 0 1], H)
    for (X, Y) in ((S, L), (L, S))
        @test !isdisjoint(X, Y)
        res, w = isdisjoint(X, Y, true)
        @test !res && w isa Vector{N} && w ∈ X && w ∈ Y
    end
    H = Hyperrectangle(N[4, 0], N[1, 1])
    L = LinearMap(N[1 0; 0 1], H)
    for (X, Y) in ((S, L), (L, S))
        @test isdisjoint(X, Y)
        res, w = isdisjoint(X, Y, true)
        @test res && w isa Vector{N} && isempty(w)
    end
end
