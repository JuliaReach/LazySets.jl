import IntervalArithmetic

# default constructor
x = Interval{Float64, IntervalArithmetic.Interval{Float64}}(IntervalArithmetic.Interval(0.0, 1.0))

# type-less constructor
x = Interval(0.0, 1.0)

@test dim(x) == 1
@test center(x) == [0.5]
@test low(x) == 0.0 && high(x) == 1.0
v = vertices_list(x)
@test [0.0] in v && [1.0] in v
# test interface method an_element and membership
@test an_element(x) ∈ x
# test containment
@test (x ⊆ x) && !(x ⊆ 0.2 * x) && (x ⊆ 2. * x)
@test issubset(x, Interval(0.0, 2.0))
@test !issubset(x, Interval(-1.0, 0.5))


# + operator (= concrete Minkowski sum of intervals)
y = Interval(-2.0, 0.5)
m = x + y
@test dim(m) == 1
@test σ([1.0], m) == [1.5]
@test σ([-1.0], m) == [-2.0]
@test low(m) == -2.0 && high(m) == 1.5
v = vertices_list(m)
@test [1.5] in v && [-2.0] in v

# difference
d = x - y
@test dim(d) == 1
@test σ([1.0], d) == [3.0]
@test σ([-1.0], d) == [-0.5]
@test low(d) == -0.5 && high(d) == 3.0
v = vertices_list(d)
@test [-0.5] in v && [3.0] in v

# product of intervals: use the * operator
p = x * y
@test dim(p) == 1
@test σ([1.0], p) == [0.5]
@test σ([-1.0], p) == [-2.0]
v = vertices_list(p)
@test [0.5] in v && [-2.0] in v

# test different arithmetic operations
r = (x + y) - (d + p)
@test low(r) == -5.5 && high(r) == 4.0

# Minkowski sum (test that we get the same results as the concrete operation)
m = x ⊕ y
@test m isa MinkowskiSum
@test dim(m) == 1
@test σ([1.0], m) == [1.5]
@test σ([-1.0], m) == [-2.0]

# cartesian product
cp = x × y
@test cp isa CartesianProduct
@test dim(cp) == 2

# conversion to hyperrectangle
h = convert(Hyperrectangle, x)
@test h isa Hyperrectangle && center(h) == radius_hyperrectangle(h) == [0.5]
