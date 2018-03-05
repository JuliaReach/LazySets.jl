using IntervalArithmetic

# default constructor
x = LazySets.Interval{Float64, IntervalArithmetic.Interval{Float64}}(IntervalArithmetic.Interval(0.0, 1.0))

# type-less constructor
x = LazySets.Interval(0.0, 1.0)

# constructor with two numbers
x = LazySets.Interval(0.0, 1.0)

@test dim(x) == 1
@test center(x) == 0.5

# Minkowski sum of intervals = usual intervals sum
y = LazySets.Interval(-2.0, 0.5)
m = x + y
@test dim(m) == 1
@test σ([1.0], m) == [1.5]
@test σ([-1.0], m) == [-2.0]

# Product of intervals: uses the * operator
p = x * y
@test dim(p) == 1
@test σ([1.0], p) == [0.5]
@test σ([-1.0], p) == [-2.0]
