function tridiagm(a, b, c, n)  # to generate tridiagonal matrices
       dd, du = ones(n), ones(n - 1)
       b*diagm(dd) + a*diagm(du, -1) + c*diagm(du, 1)
end

SUITE["Balls"] = BenchmarkGroup()
SUITE["Linear map"] = BenchmarkGroup()
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
