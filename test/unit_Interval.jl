import IntervalArithmetic

for N in [Float64, Float32, Rational{Int}]
    # constructor from IntervalArithmetic.Interval
    x = Interval(IntervalArithmetic.Interval(N(0), N(1)))

    # constructor from a vector
    x = Interval(N[0, 1])

    # type-less constructor
    x = Interval(N(0), N(1))

    @test dim(x) == 1
    @test center(x) == N[0.5]
    @test low(x) == N(0) && high(x) == N(1)
    v = vertices_list(x)
    @test N[0] in v && N[1] in v
    # test interface method an_element and membership
    @test an_element(x) ∈ x
    # test containment
    @test (x ⊆ x) && !(x ⊆ N(0.2) * x) && (x ⊆ N(2) * x)
    @test ⊆(x, Interval(N(0), N(2)))
    @test !⊆(x, Interval(N(-1), N(0.5)))

    # radius_hyperrectangle
    @test radius_hyperrectangle(x) == [N(0.5)]
    @test radius_hyperrectangle(x, 1) == N(0.5)

    # + operator (= concrete Minkowski sum of intervals)
    y = Interval(N(-2), N(0.5))
    m = x + y
    @test dim(m) == 1
    @test σ(N[1], m) == N[1.5]
    @test σ(N[-1], m) == N[-2]
    @test low(m) == N(-2) && high(m) == N(1.5)
    v = vertices_list(m)
    @test N[1.5] in v && N[-2] in v

    # difference
    d = x - y
    @test dim(d) == 1
    @test σ(N[1], d) == N[3]
    @test σ(N[-1], d) == N[-0.5]
    @test low(d) == N(-0.5) && high(d) == N(3)
    v = vertices_list(d)
    @test N[-0.5] in v && N[3] in v

    # product of intervals: use the * operator
    p = x * y
    @test dim(p) == 1
    @test σ(N[1], p) == N[0.5]
    @test σ(N[-1], p) == N[-2]
    v = vertices_list(p)
    @test N[0.5] in v && N[-2] in v

    # test different arithmetic operations
    r = (x + y) - (d + p)
    @test low(r) == N(-5.5) && high(r) == N(4)

    # isempty
    @test !isempty(x)

    # Minkowski sum (test that we get the same results as the concrete operation)
    m = x ⊕ y
    @test m isa MinkowskiSum
    @test dim(m) == 1
    @test σ(N[1], m) == N[1.5]
    @test σ(N[-1], m) == N[-2]

    # cartesian product
    cp = x × y
    @test cp isa CartesianProduct
    @test dim(cp) == 2

    # conversion to hyperrectangle
    h = convert(Hyperrectangle, x)
    @test h isa Hyperrectangle && center(h) == radius_hyperrectangle(h) == N[0.5]
end
