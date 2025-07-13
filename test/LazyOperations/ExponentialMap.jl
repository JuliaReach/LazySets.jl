for exp_backend in [ExponentialUtilities, Expokit]
    # test changing the exponentiation backend
    LazySets.set_exponential_backend!(exp_backend)

    for N in @tN([Float64, Float32])
        # dimension (choose a multiple of 3)
        n = 3 * 2

        # sparse matrix
        m = spzeros(N, n, n)
        m[1, 1] = N(1 // 25)
        m[1, 2] = N(1 // 5)
        m[3, 2] = N(-8 // 5)
        m[4, 2] = N(-3 // 2)
        m[5, 2] = N(-6 // 5)
        m[6, 3] = N(-3 // 10)
        m[2, 4] = N(1 // 2)
        m[5, 4] = N(-3 // 4)
        m[6, 4] = N(-1 // 2)
        m[3, 5] = N(-3 // 4)
        m[2, 6] = N(4 / 5)
        m[3, 6] = N(7 // 10)
        m[4, 6] = N(-1 // 2)
        m[5, 6] = N(3 // 10)

        # a set
        b = BallInf(ones(N, n), N(0.1))

        # linear map
        lm = m * b

        # the (lazy) exponential of a sparse matrix
        me = SparseMatrixExp(m)

        # constructor from a dense matrix
        @test_throws ArgumentError SparseMatrixExp(Matrix{N}(I, 2, 2))

        # constructor from a non-square matrix
        # we would use sprandn(N, 2, 3, 0.1) but it is available from Julia v1.1
        mr = randn(N, 2, 3)
        @test_throws AssertionError SparseMatrixExp(sparse(mr))

        # size
        @test size(me, 1) == n
        @test size(me) == (n, n)

        # columns & rows
        me2 = SparseMatrixExp(sparse(N(1) * I, n, n))
        v = zeros(N, n)
        v[1] = N(ℯ)
        @test get_columns(me2, [1, 2])[:, 1] == get_column(me2, 1) == v
        @test get_rows(me2, [1, 2])[1, :] == transpose(get_row(me2, 1)) == v

        # the exponential map of a convex set, same as ExponentialMap(me, b)
        emap = me * b
        @test emap isa ExponentialMap{N}
        # passing a sparse matrix automatically converts to SparseMatrixExp
        @test emap == ExponentialMap(m, b)

        # absorbing elements
        for neutral in [ZeroSet{N}(6), EmptySet{N}(6)]
            @test me * neutral == neutral
        end

        # dimension
        @test dim(emap) == n

        # the support vector of an exponential map
        d = N[2, -1, -1, 0, 1, -1]
        svec = σ(d, emap)
        # check that it works with sparse vectors
        σ(sparsevec(d), emap)

        # check consistency with respect to explicit computation of the matrix exponential
        svec_explicit = σ(d, exp(Matrix(m)) * b)
        @test svec ≈ svec_explicit

        # concretize
        X = concretize(emap)
        @test σ(d, X) ≈ svec

        # boundedness
        @test isbounded(emap) && isboundedtype(typeof(emap))
        emap2 = me * HalfSpace(ones(N, n), N(1))
        @test !isbounded(emap2) && !isboundedtype(typeof(emap2))

        # ispolyhedral
        @test ispolyhedral(emap)
        if N isa AbstractFloat
            emap2 = me * Ball2(zeros(N, n), N(1))
            @test !ispolyhedral(emap2)
        end

        # isempty
        @test !isempty(emap)

        # membership
        x = ones(N, n)
        @test x ∈ b
        @test me2 * x ∈ me2 * b
        @test -x ∉ b
        @test me2 * -x ∉ me2 * b

        # construct an exponential map where we only pass the matrix
        #ExponentialMap(m, b) # ExponentialMap

        # the exponential map of an exponential map falls back to a Cartesian product
        cpem = emap * emap # CartesianProduct of ExponentialMap

        # the exponential map applied to an exponential map
        ExponentialMap(SparseMatrixExp(2 * m), emap)

        # symmetric_interval_hull
        sih = symmetric_interval_hull(emap)
        @test sih isa Hyperrectangle{N} && center(sih) == zeros(N, 6)
        @test radius_hyperrectangle(sih) ≈
              N[1.4615602805461578, 1.7892495373142225, 1.1215454866370536,
                1.9033524001317403, 0.5475680039922208, 1.3374631184550847]

        ### ExponentialProjectionMap

        # for projection tests, assume that n is divisible by three
        @assert mod(n, 3) == 0
        nb = div(n, 3)
        # the projection of exp(A) on the (m, m)-dimensional right-most upper block
        R = [spzeros(N, nb, nb); spzeros(N, nb, nb); sparse(N(1) * I, nb, nb)]
        L = [sparse(N(1) * I, nb, nb) spzeros(N, nb, nb) spzeros(N, nb, nb)]
        proj = ProjectionSparseMatrixExp(L, me, R)

        # build an exponential projection map : it is the application of the projection
        # of exp(A) over a given set
        b = BallInf(ones(N, nb), N(1))
        projmap = proj * b

        # query the ambient dimension of projmap (hint: it is the output dimension)
        @test dim(projmap) == nb

        # boundedness
        @test isbounded(projmap)
        @test isbounded(ProjectionSparseMatrixExp(spzeros(N, nb, nb), me, R) *
                        HalfSpace(ones(N, nb), N(1)))
        @test isbounded(ProjectionSparseMatrixExp(L, me, spzeros(N, nb, nb)) *
                        HalfSpace(ones(N, nb), N(1)))
        @test !isbounded(proj * HalfSpace(ones(N, nb), N(1)))

        # isempty
        @test !isempty(projmap)

        # support vector and support function
        d = N[3 // 20, 11 // 50]
        for d in (d, sparsevec(d))  # check that it works with sparse vectors
            @test σ(d, projmap) ≈ N[0.1418094366, 1.2940902111]
            @test ρ(d, projmap) ≈ N(0.3059712619)
        end

        # check consistency with respect to explicit computation of the matrix exponential
        P = L * exp(Matrix(m)) * R
        @test σ(d, projmap) ≈ σ(d, P * b)
        @test ρ(d, projmap) ≈ ρ(d, P * b)

        # vertices_list
        b = BallInf(N[0, 0], N(1))
        M = SparseMatrixExp(spzeros(N, 2, 2))
        vlist = vertices_list(ExponentialMap(M, b))
        @test ispermutation(vlist, [N[1, 1], N[-1, 1], N[1, -1], N[-1, -1]])

        # concretize
        X = concretize(projmap)
        @test X ≈ Zonotope(N[0.070904718, 0.6470451],
                           N[0.00139082656 0.069513891; 0.026896709 0.6201484])
    end
end
