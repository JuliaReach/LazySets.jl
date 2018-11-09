for N in [Float64, Rational{Int}, Float32]
    # random line
    rand(Line)

    # construction
    l1 = Line(N[0, 1], N(1))
    l2 = Line(N[1, 0], N(1))
    l3 = Line(N[0, 1], N(2))
    l4 = Line(N[1, 1], N(0))
    l1_copy = Line(N[0, 1], N(1))

    # alternative construction from two points
    # line with positive slope
    L = Line([1.0, 0.0], [2.0, 0.5])
    @test [1.0, 0.0] ∈ L &&
          [2.0, 0.5] ∈ L &&
          [2.0, 1.5] ∉ L

    # vertical line
    L = Line([1.0, 0.0], [1.0, 2.5])
    @test [1.0, 0.0] ∈ L &&
          [1.0, 2.5] ∈ L &&
          [1.0, 10.0] ∈ L
          [10.0, 10.0] ∉ L

    # dimension
    @test dim(l1) == 2

    # support vector
    σ(N[0, 1], l1)
    σ(N[1, 0], l2)
    σ(N[0, 1], l3)

    # isempty
    @test !isempty(l1)

    # an_element and membership
    an_element(l1) ∈ l1
    an_element(l2) ∈ l2
    an_element(l3) ∈ l3
    an_element(l4) ∈ l4

    # constrained dimensions
    @test constrained_dimensions(l1) == [2]
    @test constrained_dimensions(l4) == [1, 2]

    # concrete intersection
    cap11 = intersection(l1, l1_copy)
    cap12 = intersection(l1, l2)
    cap13 = intersection(l1, l3)
    @test cap11 isa Line && cap11.a == l1.a && cap11.b == l1.b
    @test cap12 isa Singleton && element(cap12) == N[1, 1]
    @test cap13 isa EmptySet{N}
end
