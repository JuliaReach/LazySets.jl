using LazySets.Arrays: _rationalize

for N in [Float32, Float64]
    x = [N(1), N(2), N(3)]
    out = [1//1, 2//1, 3//1]
    @test _rationalize(x) == out
    v = _rationalize(BigInt, x)
    @test v[1] isa Rational{BigInt}
    @test _rationalize(x, tol=2*eps(N)) == out
end
