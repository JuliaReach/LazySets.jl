using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # compact set
    @test Ball1 <: CompactSet && !(Ball1 <: NonCompactSet)
    b = Ball1(N[0, 1], N(1))
    @test b isa CompactSet && !(b isa NonCompactSet)

    # non-compact set
    @test Line2D <: NonCompactSet && !(Line2D <: CompactSet)
    l = Line2D(N[0, 1], N(1))
    @test l isa NonCompactSet && !(l isa CompactSet)

    # lazy operations are neither compact nor non-compact by default
    @test !(MinkowskiSum <: CompactSet) && !(MinkowskiSum <: NonCompactSet)
    ms = b + l
    @test !(ms isa CompactSet) && !(ms isa NonCompactSet)
end
