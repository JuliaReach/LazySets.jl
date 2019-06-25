# generate a tridiagonal matrix
function tridiagm(a, b, c, n)
    dd, du = ones(n), ones(n - 1)
    b*diagm(0 => dd) + a*diagm(-1 => du) + c*diagm(1 => du)
end

# support function
SUITE["ρ"] = BenchmarkGroup()

# support vector
SUITE["σ"] = BenchmarkGroup()

N = Float64
for set_type in (Ball1, BallInf, Hyperrectangle)
    for n in (2, 10, 100, 1000)
        rng = MersenneTwister(n)
        d = rand(rng, n)
        X = rand(set_type, dim=n, rng=rng)
        SUITE["ρ"][string(set_type), "dense", n] = @benchmarkable ρ($d, $X)
        SUITE["σ"][string(set_type), "dense", n] = @benchmarkable σ($d, $X)
        A = tridiagm(1, -2, 0.5, n)
        Y = A * X
        SUITE["ρ"]["LinearMap of $(string(set_type))", "dense", n] = @benchmarkable ρ($d, $Y)
        SUITE["σ"]["LinearMap of $(string(set_type))", "dense", n] = @benchmarkable σ($d, $Y)
    end
end
