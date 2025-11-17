using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    E = Intersection(HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-1)))  # empty set
    Xnc = UnionSet(BallInf(N[1], N(1)), BallInf(N[4], N(1)))  # nonconvex set

    # linear_combination, 1D case with empty set
    for Y in (linear_combination(E, Xnc), linear_combination(Xnc, E))
        @test Y == E
    end
end
