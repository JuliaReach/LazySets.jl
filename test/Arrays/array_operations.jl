using LazySets.Arrays: _rationalize, append_zeros, prepend_zeros

for N in [Float64, Float32, Rational{Int}]
    v1 = N[0, 4, 0]
    v2 = sparsevec(N[2], N[4], 3)
    v3 = SingleEntryVector(2, 3, N(4))
    @test v1 == v2 == v3
    v1a = append_zeros(v1, 2)
    v1p = prepend_zeros(v1, 2)
    @test v1a == append_zeros(v2, 2) == append_zeros(v3, 2) == N[0, 4, 0, 0, 0]
    @test v1p == prepend_zeros(v2, 2) == prepend_zeros(v3, 2) == N[0, 0, 0, 4, 0]
end

for N in [Float32, Float64]
    x = [N(1), N(2), N(3)]
    out = [1//1, 2//1, 3//1]
    @test _rationalize(x) == out
    v = _rationalize(BigInt, x)
    @test v[1] isa Rational{BigInt}
    @test _rationalize(x, tol=2*eps(N)) == out
end
