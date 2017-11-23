∅ = EmptySet()
B = BallInf(ones(2), 0.5)

# absorbing for the Cartesian product
@test B * ∅ isa EmptySet && ∅ * B isa EmptySet

# absorbing for the Minkowski sum
@test B + ∅ isa EmptySet && ∅ + B isa EmptySet
