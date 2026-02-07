using LazySets, Test
IA = LazySets.IA
using LazySets.IA: IntervalBox
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

for N in [Float64]
    # rounding error
    P = VPolytope([[-0.05743585962166445, 0.1076870835851457], [-0.05449856438167133, 0.09267590692998105], [-0.0469929760540887, 0.08810761384239148], [-0.046904863337052544, 0.08811328836520504], [-0.046904809759468036, 0.08811334365928523], [-0.0498270983945782, 0.10311537329191771], [-0.04982715636089584, 0.1031155435919711], [-0.05732817117721204, 0.10769426977317734]])
    H = Hyperrectangle([-0.5, 0.5], [0.1, 0.1])
    @test isdisjoint(P, H)
end
