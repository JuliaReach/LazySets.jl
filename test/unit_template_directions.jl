for N in [Float64, Float32, Rational{Int}]
    # unit vector
    for n in 1:3
        for i in 1:n
            e = SingleEntryVector(i, n, one(N))
            for j in 1:n
                @test e[j] == (i == j ? one(N) : zero(N))
            end
        end
    end

    # check that #1359 is fixed
    e = SingleEntryVector(2, 3, N(1)) # vector [0, 1, 0]
    D = Diagonal(N[1, 2, 3])
    De = D * e
    @test De isa SingleEntryVector && Vector(De) == N[0, 2, 0]

    # template direction approximation
    for n in 1:3
        B = BallInf(zeros(N, n), N(2))
        A = Matrix{N}(2I, n, n) + ones(N, n, n)
        X = A * B

        # box directions
        dir = BoxDirections{N}(n)
        @test isbounding(dir)
        @test dim(dir) == n
        box = overapproximate(X, dir)
        @test length(dir) == length(box.constraints) == 2*n
        box = overapproximate(X, BoxDirections)
        @test length(dir) == length(box.constraints) == 2*n

        # octagon directions
        dir = OctDirections{N}(n)
        @test isbounding(dir)
        @test dim(dir) == n
        oct = overapproximate(X, dir)
        @test length(dir) == length(oct.constraints) == 2 * n^2
        oct = overapproximate(X, OctDirections)
        @test length(dir) == length(oct.constraints) == 2 * n^2

        # box-diagonal directions
        dir = BoxDiagDirections{N}(n)
        @test isbounding(dir)
        @test dim(dir) == n
        boxdiag = overapproximate(X, dir)
        @test length(dir) == length(boxdiag.constraints) ==
              (n == 1 ? 2 : 2^n + 2*n)
        boxdiag = overapproximate(X, BoxDiagDirections)
        @test length(dir) == length(boxdiag.constraints) ==
              (n == 1 ? 2 : 2^n + 2*n)

        # spherical directions approximation
        if n == 2 && N in [Float32, Float64]
            dir = PolarDirections{N}(2)
            @test !isbounding(dir)
            dir = PolarDirections{N}(5)
            @test isbounding(dir)
            @test dim(dir) == 2
            polar = overapproximate(X, dir)
        end

        # spherical directions approximation
        if n == 3 && N in [Float32, Float64]
            dir = SphericalDirections{N}(2, 2)
            @test !isbounding(dir)
            dir = SphericalDirections{N}(5, 5)
            @test isbounding(dir)
            @test dim(dir) == 3
            spherical = overapproximate(X, dir)
        end

        # overapproximate lazy polyhedral intersections
        if N in [Float64]
            Y = B ∩ Ball1(zeros(N, n), N(1))
            Z = B ∩ HalfSpace(ones(N, n), N(1))
            for dir in [BoxDirections, OctDirections, BoxDiagDirections]
                overapproximate(Y, dir)
                overapproximate(Z, dir)
            end
        end
    end
end

for N in [Float64]
    for n in 1:3
        B = BallInf(zeros(N, n), N(2))
        A = Matrix{N}(2I, n, n) + ones(N, n, n)
        X = A * B

        # custom directions
        # empty list of directions
        dir = CustomDirections(Vector{N}[]; n=n)
        @test dim(dir) == n
        P = overapproximate(X, dir)
        @test isempty(constraints_list(P))
        # minimal number of bounded directions
        if n == 1
            dirs = [N[-1], N[1]]
        elseif n == 2
            dirs = [N[-1, 0], N[0, -1], N[1, 1]]
        elseif n == 3
            dirs = [N[-1, 0, 0], N[0, -1, 0], N[0, 0, -1], N[1, 1, 1]]
        end
        dir = CustomDirections(dirs)
        @test isbounding(dir)
        @test !isbounding(CustomDirections(dirs[1:end-1]))
        P = overapproximate(X, dir)
        @test P isa HPolytope && length(constraints_list(P)) == length(dirs)
    end
end

# default Float64 constructors
BoxDirections(3)
OctDirections(3)
BoxDiagDirections(3)
SphericalDirections(3)
