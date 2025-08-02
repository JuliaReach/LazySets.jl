using LazySets, Test
using LazySets.ReachabilityBase.Arrays: ispermutation
IA = LazySets.IA
using LazySets.IA: IntervalBox

for N in @tN([Float64, Float32, Rational{Int}])
    # vertices_list for IA types
    Y = IntervalBox(IA.interval(N(0), N(1)), IA.interval(N(2), N(3)))
    res = vertices_list(Y)
    @test ispermutation(res, [[N(1), N(3)], [N(0), N(3)], [N(1), N(2)], [N(0), N(2)]])
    res = vertices_list(Y[1])
    @test ispermutation(res, [[N(0)], [N(1)]])
end
