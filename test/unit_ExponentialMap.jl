for N in [Float64, Float32]
    # dimension (choose a multiple of 3)
    n = 3*2

#     p = 0.4 # occupation probability for these benchs
#     m = sprandn(n, n, p)

    # sparse matrix
    m = spzeros(N, n, n)
    m[1, 1] = to_N(N, 0.0439726)
    m[1, 2] = to_N(N, 0.206197)
    m[3, 2] = to_N(N, -1.64948)
    m[4, 2] = to_N(N, -1.5343)
    m[5, 2] = to_N(N, -1.22162)
    m[6, 3] = to_N(N, -0.281368)
    m[2, 4] = to_N(N, 0.535221)
    m[5, 4] = to_N(N, -0.768595)
    m[6, 4] = to_N(N, -0.499827)
    m[3, 5] = to_N(N, -0.76484)
    m[2, 6] = to_N(N, 0.821664)
    m[3, 6] = to_N(N, 0.715326)
    m[4, 6] = to_N(N, -0.545632)
    m[5, 6] = to_N(N, 0.312998)

    # a set
    b = BallInf(ones(N, n), N(0.1))

    # linear map
    lm = m * b

    # the (lazy) exponential of a sparse matrix
    me = SparseMatrixExp(m)

    # constructor from a dense matrix
    @test_throws ErrorException SparseMatrixExp(eye(N, 2))

    # size
    @test size(me, 1) == n
    @test size(me) == (n, n)
    # product of the exponential maps of two commuting matrices
    # WARNING: assuming commutativity of matrix exponents
    #me * me

    # check that the eye method is defined
    @test hasmethod(eye, Tuple{typeof(me)})
    y = to_N(N, eye(me))
    @test norm(y - diagm(ones(N, n))) == 0

    # columns & rows
    me2 = SparseMatrixExp(sparse(N(1) * I, n, n))
    v = zeros(N, n)
    v[1] = N(ℯ)
    @test get_columns(me2, [1, 2])[:, 1] == get_column(me2, 1) == v
    @test get_rows(me2, [1, 2])[1, :] == transpose(get_row(me2, 1)) == v

    # the exponential map of a convex set, same as ExponentialMap(me, b)
    emap = me * b
    @test emap isa ExponentialMap{N}

    # dimension
    @test dim(emap) == n

    # the support vector of an exponential map
#     d = randn(N, n)
    d = to_N(N, [2.29681, -0.982841, -0.642168, 0.0167593, 1.32862, -0.855418])
    svec = σ(d, emap)
    # check that it works with sparse vectors
    σ(sparsevec(d), emap)

    # check consistency with respect to explicit computation of the matrix exponential
    svec_explicit = σ(d, expm(full(m)) * b)
    if N == Float64
        # precision with Float32 is not sufficient
        @test svec ≈ svec_explicit
    end

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
    ExponentialMap(SparseMatrixExp(2*m), emap)

    # for projection tests, let's assume that n is divisible by three
    assert(mod(n, 3) == 0)
    nb = div(n, 3)
    # the projection of exp(A) on the (m, m)-dimensional right-most upper block
    R = [spzeros(N, nb, nb); spzeros(N, nb, nb); sparse(N(1) * I, nb, nb)]
    L = [sparse(N(1) * I, nb, nb) spzeros(N, nb, nb) spzeros(N, nb, nb)]
    proj = ProjectionSparseMatrixExp(L, me, R)

    # build an exponential projection map : it is the application of the projection
    # of exp(A) over a given set
    b = BallInf(ones(N, nb), to_N(N, 0.1))
    projmap = proj * b

    # query the ambient dimension of projmap (hint: it is the output dimension)
    @test dim(projmap) == nb

    #compute the support vector of the projection of an exponential map
#     d = randn(N, nb)
    d = to_N(N, [0.152811, 0.22498])
    svec = σ(d, projmap)
    # check that it works with sparse vectors
    σ(sparsevec(d), projmap)

    # check consistency with respect to explicit computation of the matrix exponential
    P = L * expm(full(m)) * R
    svec_explicit = σ(d, P*b)
    @test svec ≈ svec_explicit
end
