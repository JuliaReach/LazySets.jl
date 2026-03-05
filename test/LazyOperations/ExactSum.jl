using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    P1 = SparsePolynomialZonotope(N[4, 4], N[2 1 2; 0 2 2], hcat(N[1; 0]), [1 0 3; 0 1 1], [1, 3])
    P2 = SparsePolynomialZonotope(zeros(N, 2), N[2 0 1; 1 2 1], zeros(N, 2, 0), [1 0 1; 0 1 3], [1, 3])

    ES = ExactSum(P1, P2)
    @test ES.X == P1
    @test ES.Y == P2

    @test P1 ⊞ P2 == ES
    @test isoperationtype(typeof(ES))
    @test LazySets.concrete_function(ExactSum) == exact_sum
    @test concretize(ES) == exact_sum(P1, P2)
end
