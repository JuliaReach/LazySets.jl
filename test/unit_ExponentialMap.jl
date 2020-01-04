for N in [Float64, Float32]
    # dimension (choose a multiple of 3)
    n = 3*2

    # sparse matrix
    m = spzeros(N, n, n)
    m[1, 1] = N(1//25)
    m[1, 2] = N(1//5)
    m[3, 2] = N(-8//5)
    m[4, 2] = N(-3//2)
    m[5, 2] = N(-6//5)
    m[6, 3] = N(-3//10)
    m[2, 4] = N(1//2)
    m[5, 4] = N(-3//4)
    m[6, 4] = N(-1//2)
    m[3, 5] = N(-3//4)
    m[2, 6] = N(4/5)
    m[3, 6] = N(7//10)
    m[4, 6] = N(-1//2)
    m[5, 6] = N(3//10)

    # a set
    b = BallInf(ones(N, n), N(0.1))

    # linear map
    lm = m * b

    # the (lazy) exponential of a sparse matrix
    me = SparseMatrixExp(m)

    # constructor from a dense matrix
    @test_throws ErrorException SparseMatrixExp(Matrix{N}(I, 2, 2))

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

    # boundedness
    @test isbounded(emap)
    @test !isbounded(me * HalfSpace(ones(N, n), N(1)))

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
    ExponentialMap(SparseMatrixExp(2*m), emap)

    # for projection tests, let's assume that n is divisible by three
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
    @test isbounded(ProjectionSparseMatrixExp(spzeros(N, nb, nb), me, R) * HalfSpace(ones(N, nb), N(1)))
    @test isbounded(ProjectionSparseMatrixExp(L, me, spzeros(N, nb, nb)) * HalfSpace(ones(N, nb), N(1)))
    # the following test crashes because ρ(::ExponentialProjectionMap) is not implemented yet
    @test_throws ErrorException !isbounded(proj * HalfSpace(ones(N, nb), N(1)))

    # isempty
    @test !isempty(projmap)

    #compute the support vector of the projection of an exponential map
    d = N[3//20, 11//50]
    svec = σ(d, projmap)
    # check that it works with sparse vectors
    svec == σ(sparsevec(d), projmap)

    # check consistency with respect to explicit computation of the matrix exponential
    P = L * exp(Matrix(m)) * R
    svec_explicit = σ(d, P*b)
    @test svec ≈ svec_explicit

    # vertices_list
    b = BallInf(N[0, 0], N(1))
    M = SparseMatrixExp(spzeros(N, 2, 2))
    vlist = vertices_list(ExponentialMap(M, b))
    @test ispermutation(vlist, [N[1, 1], N[-1, 1], N[1, -1], N[-1, -1]])
end
