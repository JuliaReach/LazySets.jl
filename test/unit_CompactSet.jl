for N in [Float64, Rational{Int}, Float32]
    # compact set
    @test Ball1 <: CompactSet && !(Ball1 <: NonCompactSet)
    b = Ball1(N[0, 1], N(1))
    @test b isa CompactSet && !(b isa NonCompactSet)

    # non-compact set
    @test Line <: NonCompactSet && !(Line <: CompactSet)
    l = Line(N[0, 1], N(1))
    @test l isa NonCompactSet && !(l isa CompactSet)

    # lazy operations are neither compact nor non-compact by default
    @test !(MinkowskiSum <: CompactSet) && !(MinkowskiSum <: NonCompactSet)
    ms = b + l
    @test !(ms isa CompactSet) && !(ms isa NonCompactSet)
end
