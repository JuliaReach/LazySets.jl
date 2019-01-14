for N in [Float64, Rational{Int}, Float32]
    B1 = BallInf(zeros(N, 2), N(1))
    B2 = Ball1(ones(N, 2), N(1))
    UXY = UnionSet(B1, B2)

    # type alias
    U = B1 ∪ B2

    # array type (union of a finite number of convex sets)
    Uarr = UnionSetArray([B1, B2])

    for U in [UXY, Uarr]
        # support vector (default algorithm)
        d = N[1, 0]
        @test σ(d, U) == [N(2), N(1)]
        @test σ(d, U, algorithm="support_vector") == [N(2), N(1)]

        # support vector (support function algorithm)
        @test σ(d, U, algorithm="support_function") == [N(2), N(1)]

        # support function
        @test ρ(d, U) == N(2)

        # an_element and membership
        @test an_element(U) ∈ U

        # emptiness
        @test !isempty(U)
    end

    # emptiness
    emptyP = HPolyhedron([HalfSpace(N[1, 0], N(0)), HalfSpace(N[-1, 0], N(-1))])
    @test !isempty(emptyP ∪ B1) && !isempty(B1 ∪ emptyP) &&
          isempty(emptyP ∪ emptyP)
    @test !isempty(UnionSetArray([emptyP, B1])) &&
          !isempty(UnionSetArray([B1, B2, emptyP])) &&
          isempty(UnionSetArray([emptyP, emptyP]))
end
