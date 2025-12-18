using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    H1 = BallInf(N[0, 0], N(1))
    H2 = BallInf(N[0, 0], N(2))

    # strict inclusion
    for (X, Y) in [(H1, H2)]
        res, w = ⊂(X, Y, true)
        @test X ⊂ Y && res && w ∈ Y && w ∉ X
    end
    # no strict inclusion
    res, w = ⊂(H2, H1, true)
    @test !(H2 ⊂ H1) && !res && w isa Vector{N} && w ∈ H2 && w ∉ H1
    # equal sets
    res, w = ⊂(H1, H1, true)
    @test !(H1 ⊂ H1) && !res && w isa Vector{N} && isempty(w)
end
