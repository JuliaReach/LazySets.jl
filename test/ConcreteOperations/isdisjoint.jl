using LazySets, Test
IA = LazySets.IA
using LazySets.IA: IntervalBox
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
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
end
