using LazySets: translate, ispermutation

for N in [Float64, Rational{Int}, Float32]
    # random line
    rand(Line)

    # construction
    l1 = Line(from=N[0, 1], to=N[1, 1]) # two points in the line
    l2 = Line(N[0, 1], N[1, 0]) # point and direction

    # the lines are the same modulo the sign of the normal vector
    @test l1.p ≈ l2.p && l1.n ≈ -l2.n

    # dimension
    @test dim(l1) == 2

    # boundedness
    @test !isbounded(l1)

    # universality
    @test !isuniversal(l1)

    # isempty
    @test !isempty(l1)

    # an_element and membership
    an_element(l1) ∈ l1

    # constrained dimensions
    @test constrained_dimensions(l1) == [2]

    # constraints_list
    if N <: AbstractFloat
        @test ispermutation(constraints_list(l2),
                            [HalfSpace(N[0, 1], N(1)),    # y <= 1
                             HalfSpace(N[0, -1], N(-1))]) # y >= 1
    end

    # translation
    @test translate(l2, N[0, 1]) == Line(N[0, 2], N[1, 0])
end
