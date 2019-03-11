for N in [Float64, Rational{Int}, Float32]
    B1 = BallInf(zeros(N, 2), N(1))
    B2 = Ball1(ones(N, 2), N(1))
    B3 = Hyperrectangle(low=N[-1, -1], high=N[2, 2])
    S = Singleton(N[10, 10])
    UXY = UnionSet(B1, B2)

    # type alias
    @test B1 ∪ B2 == UXY

    # array type (union of a finite number of convex sets)
    Uarr = UnionSetArray([B1, B2])

    for U in [UXY, Uarr]
        # dimension
        @test dim(U) == dim(B1)

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

        # boundedness
        @test isbounded(U)

        # inclusion
        subset, point = ⊆(U, B3, true)
        @test U ⊆ B3 && subset && point == N[]
        subset, point = ⊆(U, B2, true)
        @test !(U ⊆ B2) && !subset && point ∈ U && point ∉ B2

        # isdisjoint
        disjoint1, point1 = isdisjoint(U, B1, true)
        disjoint2, point2 = isdisjoint(B1, U, true)
        @test !isdisjoint(U, B1) && !isdisjoint(B1, U) && !disjoint1 &&
              !disjoint2 && point1 ∈ U && point1 ∈ B1 && point2 ∈ U &&
              point2 ∈ B1
        disjoint1, point1 = isdisjoint(U, S, true)
        disjoint2, point2 = isdisjoint(S, U, true)
        @test isdisjoint(U, S) && isdisjoint(S, U) && disjoint1 &&
              disjoint2 && point1 == point2 == N[]
    end

    # emptiness
    emptyP = HPolyhedron([HalfSpace(N[1, 0], N(0)), HalfSpace(N[-1, 0], N(-1))])
    @test !isempty(emptyP ∪ B1) && !isempty(B1 ∪ emptyP) &&
          isempty(emptyP ∪ emptyP)
    @test !isempty(UnionSetArray([emptyP, B1])) &&
          !isempty(UnionSetArray([B1, B2, emptyP])) &&
          isempty(UnionSetArray([emptyP, emptyP]))

    # boundedness
    unboundedP = HPolyhedron([HalfSpace(N[1, 0], N(0))])
    @test !isbounded(unboundedP ∪ B1) && !isbounded(B1 ∪ unboundedP)
    @test !isbounded(UnionSetArray([unboundedP, B1])) &&
          !isbounded(UnionSetArray([B1, B2, unboundedP]))

    # tests that only work with Float64
    if N in [Float64]
        for U in [UXY, Uarr]
            # intersection
            @test !isempty(intersection(U, B3)) && !isempty(intersection(B3, U))
            @test isempty(intersection(U, S)) && isempty(intersection(S, U))
        end

        # isdisjoint
        disjoint1, point1 = isdisjoint(UXY, UXY, true)
        disjoint2, point2 = isdisjoint(UXY, Uarr, true)
        disjoint3, point3 = isdisjoint(Uarr, UXY, true)
        disjoint4, point4 = isdisjoint(Uarr, Uarr, true)
        @test !isdisjoint(UXY, UXY) && !isdisjoint(UXY, Uarr) &&
              !isdisjoint(Uarr, UXY) && !isdisjoint(Uarr, Uarr) && !disjoint1 &&
              !disjoint2 && !disjoint3 && !disjoint4 && point1 ∈ UXY &&
              point2 ∈ UXY && point2 ∈ Uarr && point3 ∈ UXY && point3 ∈ Uarr &&
              point4 ∈ Uarr
    end
end
