import LazySets.Approximations,
       Approximations.UnitVector,
       Approximations.BoxDirections,
       Approximations.OctDirections,
       Approximations.BoxDiagDirections,
       Approximations.overapproximate

for N in [Float64, Float32, Rational{Int}]
    # unit vector
    for n in 1:3
        for i in 1:n
            e = UnitVector(i, n, one(N))
            for j in 1:n
                @test e[j] == (i == j ? one(N) : zero(N))
            end
        end
    end

    # template direction approximation
    for n in 1:3
        B = BallInf(zeros(N, n), N(2.))
        A = 2.*eye(N, n) + ones(N, n, n)
        X = A * B

        # box directions
        dir = BoxDirections{N}(n)
        @test dim(dir) == n
        box = overapproximate(X, dir)
        @test length(dir) == length(box.constraints) == 2*n

        # octagon directions
        dir = OctDirections{N}(n)
        @test dim(dir) == n
        oct = overapproximate(X, dir)
        @test length(dir) == length(oct.constraints) == 2 * n^2

        # box-diagonal directions
        dir = BoxDiagDirections{N}(n)
        @test dim(dir) == n
        boxdiag = overapproximate(X, dir)
        @test length(dir) == length(boxdiag.constraints) ==
              (n == 1 ? 2 : 2^n + 2*n)
    end
end

# default Float64 constructors
BoxDirections(3)
OctDirections(3)
BoxDiagDirections(3)
