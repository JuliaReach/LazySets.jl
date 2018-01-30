for N in [Float64] # TODO Float32
    # dimension (choose a multiple of 3)
    n = 3*2

    # occupation probability for these benchs
    p = 0.4

    # sparse matrix
    m = sprandn(n, n, p)
    m = convert(SparseMatrixCSC{N,Int}, m)

    # a set
    b = BallInf(ones(N, n), to_N(N, 0.1))

    # linear map
    lm = m * b

    # the (lazy) exponential of a sparse matrix
    me = SparseMatrixExp(m)

    # size
    @test size(me, 1) == n
    @test size(me) == (n, n)
    # product of the exponential maps of two commuting matrices
    # WARNING: assuming commutativity of matrix exponents
    #me * me

    # the exponential map of a convex set, same as ExponentialMap(me, b)
    emap = me * b

    # the support vector of an exponential map
    d = randn(N, n)
    svec = σ(d, emap)
    # check that it works with sparse vectors
    σ(sparsevec(d), emap)

    # check consistency with respect to explicit computation of the matrix exponential
    svec_explicit = σ(d, expm(full(m)) * b)
    @test svec ≈ svec_explicit

    # construct an exponential map where we only pass the matrix
    #ExponentialMap(m, b) # ExponentialMap

    # the exponential map of an exponential map falls back to a Cartesian product
    cpem = emap * emap # CartesianProduct of ExponentialMap

    # the exponential map applied to an exponential map
    emap2 = ExponentialMap(SparseMatrixExp(2*m), emap)

    # for projection tests, let's assume that n is divisible by three
    assert(mod(n, 3) == 0)
    nb = div(n, 3)
    # the projection of exp(A) on the (m, m)-dimensional right-most upper block
    R = [spzeros(N, nb, nb); spzeros(N, nb, nb); speye(N, nb, nb)]
    L = [speye(N, nb, nb) spzeros(N, nb, nb) spzeros(N, nb, nb)]
    proj = ProjectionSparseMatrixExp(L, me, R)

    # build an exponential projection map : it is the application of the projection
    # of exp(A) over a given set
    b = BallInf(ones(N, nb), to_N(N, 0.1))
    projmap = proj * b

    # query the ambient dimension of projmap (hint: it is the output dimension)
    @test dim(projmap) == nb

    #compute the support vector of the projection of an exponential map
    d = randn(N, nb)
    svec = σ(d, projmap)
    # check that it works with sparse vectors
    σ(sparsevec(d), projmap)

    # check consistency with respect to explicit computation of the matrix exponential
    P = L * expm(full(m)) * R
    svec_explicit = σ(d, P*b)
    @test svec ≈ svec_explicit
end
