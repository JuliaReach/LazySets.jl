# generate a tridiagonal matrix
function tridiagm(a, b, c, n)
    dd, du = ones(n), ones(n - 1)
    b*diagm(0 => dd) + a*diagm(-1 => du) + c*diagm(1 => du)
end

# ===========================
# Support vectors (σ)
# ===========================
SUITE["σ"] = BenchmarkGroup()

SUITE["σ"]["Balls"] = BenchmarkGroup()
SUITE["σ"]["Linear map"] = BenchmarkGroup()
begin
    for set_type in (Ball1, Ball2, BallInf)
        for n in (1, 10, 100, 1000)
            X = set_type(ones(n), 0.1)
            d = [i/n for i in 1:n]
            SUITE["Balls"][string(set_type), n] = @benchmarkable σ($d, $X)

            A = tridiagm(1, -2, 0.5, n)
            Y = A * X
            SUITE["Linear map"][string(set_type), n] = @benchmarkable σ($d, $Y)
        end
    end
end

# ===========================
# Support functions (ρ)
# ===========================
SUITE["ρ"] = BenchmarkGroup()

#SUITE["ρ"]["Ball1"] = BenchmarkGroup()
#SUITE["ρ"]["BallInf"] = BenchmarkGroup()
#SUITE["ρ"]["Hyperrectangle"] = BenchmarkGroup()
#SUITE["ρ"]["LinearMap"] = BenchmarkGroup()

for set_type in (Ball1, BallInf, Hyperrectangle)
    for n in (2, 10, 100, 1000)
        d = [i/n for i in 1:n]
        B = Ball1(ones(n), 0.1)
        X = convert(set_type, B)
        SUITE["ρ"][string(set_type), n] = @benchmarkable ρ($d, $X)
        A = tridiagm(1, -2, 0.5, n)
        Y = A * X
        SUITE["ρ"]["LinearMap", n] = @benchmarkable ρ($d, $Y)
    end
end


end
for n in (1, 10, 100, 1000)
    d = [i/n for i in 1:n]

    X = Ball1(ones(n), 0.1)
    SUITE["ρ"]["Ball1", n] = @benchmarkable ρ($d, $X)

    X = BallInf(ones(n), 0.1)
    SUITE["ρ"]["BallInf", n] = @benchmarkable ρ($d, $X)

    X = convert(Hyperrectangle, X)
    SUITE["ρ"]["Hyperrectangle", n] = @benchmarkable ρ($d, $X)

    A = tridiagm(1, -2, 0.5, n)
    Y = A * X
    SUITE["ρ"]["LinearMap", n] = @benchmarkable ρ($d, $Y)
end
