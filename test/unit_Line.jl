for N in [Float64, Rational{Int}, Float32]
    # construction
    l1 = Line(N[0, 1], N(1))
    l2 = Line(N[1, 0], N(1))
    l3 = Line(N[0, 1], N(2))
    l4 = Line(N[1, 1], N(0))
    l1_copy = Line(N[0, 1], N(1))

    # dimension
    @test dim(l1) == 2

    # support vector & membership
    σ(N[0, 1], l1)
    σ(N[1, 0], l2)
    σ(N[0, 1], l3)

    # an_element and membership
    an_element(l1) ∈ l1
    an_element(l2) ∈ l2
    an_element(l3) ∈ l3
    an_element(l4) ∈ l4

    # concrete intersection
    cap11 = intersection(l1, l1_copy)
    cap12 = intersection(l1, l2)
    cap13 = intersection(l1, l3)
    @test cap11 isa Line && cap11.a == l1.a && cap11.b == l1.b
    @test cap12 isa Singleton && element(cap12) == N[1, 1]
    @test cap13 isa EmptySet{N}
end
