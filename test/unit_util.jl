using LinearAlgebra: Diagonal
for _dummy_ in 1:1 # avoid global variable warnings
    # reseeding with random seed
    rng = LazySets.GLOBAL_RNG
    seed = rand(1:10000)
    LazySets.reseed(rng, seed)
    n1 = rand(Int)
    LazySets.reseed(rng, seed)
    n2 = rand(Int)
    @test n1 == n2

    # StrictlyIncreasingIndices
    vectors = Vector{AbstractVector{Int}}()
    for v in LazySets.StrictlyIncreasingIndices(5, 4)
        push!(vectors, copy(v))
    end
    @test vectors == [[1, 2, 3, 4], [1, 2, 3, 5], [1, 2, 4, 5], [1, 3, 4, 5], [2, 3, 4, 5]]

    # square matrix
    @test LazySets.issquare([2 3; 0 0])
    @test LazySets.issquare(sparse([1], [1], [1], 3, 3))
    @test !LazySets.issquare(hcat([2 3]))
    @test !LazySets.issquare(sparse([1], [1], [1], 3, 4))

    # invertible matrix: (1) invertible, (2) non-invertible, (3) non-square
    # dense matrix
    @test LazySets.isinvertible([2 3; 1 2])
    @test !LazySets.isinvertible([2 3; 0 0])
    @test !LazySets.isinvertible(hcat([2 3]))
    # sparse matrix
    @test LazySets.isinvertible(sparse([2 3; 1 2]))
    @test !LazySets.isinvertible(sparse([2 3; 0 0]))
    @test !LazySets.isinvertible(sparse(hcat([2 3])))
    # diagonal matrix
    @test LazySets.isinvertible(Diagonal([2 0; 0 2]))
    @test !LazySets.isinvertible(Diagonal([2 0; 0 0]))

    # matrix rank
    A = sprandn(2, 10, 0.4)
    @test rank(Matrix(A)) == rank(A) == rank(view(A, :, :))

    for N in [Float64, Rational{Int}, Float32]
        # substitution
        x = N[1, 2, 3]
        substitution = Dict(1 => N(4), 3 => N(0))
        @test LazySets.substitute(substitution, x) == N[4, 2, 0]
        LazySets.substitute!(substitution, x)
        @test x == N[4, 2, 0]

        A = N[1 4; 2 5; 3 6]
        x1 = N[0, 2, 0]
        y1 = N[3, 0]
        x2 = SingleEntryVector(2, 3, N(2))
        y2 = SingleEntryVector(1, 2, N(3))
        @test -x2 == SingleEntryVector(2, 3, N(-2))
        @test inner(x1, A, y1) == dot(x1, A * y1) == inner(x2, A, y2) ==
              dot(x2, A * y2) == N(12)

        x = N[0, 1, -1]
        @test LazySets.rectify(x) == N[0, 1, 0]

        # approximate permutation check
        v1 = [N[1, 2], N[2, 3], N[3, 4]]
        tol = Base.rtoldefault(N)/2
        v2 = [N[3 + tol, 4 - tol], N[2 - tol, 3 + tol], N[1 + tol, 2 + tol]]
        v3 = [N[2, 3], N[1, 2], N[4, 4]]
        @test ispermutation(v1, v2) && !ispermutation(v2, v3)

        # to_negative_vector
        for v in [N[-1, 0, 1], sparsevec([1, 3], N[-1, 1], 3), N(-1):N(1):N(1)]
            u = LazySets.Arrays.to_negative_vector(v)
            @test u isa Vector{N} && u == N[1, 0, -1]
        end
    end

    for N in [Float64, Float32]
        # modified dot product
        @test isnan(dot(N[1, 0], N[Inf, -Inf]))
        @test LazySets.dot_zero(N[1, 0], N[Inf, -Inf]) == N(Inf)

        # remove zero columns
        m = 3
        n = 10
        A = rand(N, m, n)
        A[:, 1] = A[:, 5] = A[:, n] = zeros(N, m)
        B = LazySets.delete_zero_columns!(A)
        @test size(B) == (m, n-3)
    end

    # Check if sets are equal
    n = 3
    I = Diagonal(ones(n))
    X = HPolytope([I; -I], ones(2*n))
    Y = convert(VPolytope, X)
    @test isequivalent(X, Y)

end
