E = EmptySet()
B = BallInf(ones(2), 1.)

# testing that the empty set is an absorbing element for the cartesian product
@test B * E isa EmptySet && E * B isa EmptySet
# testing the mathematical alias ×
@test B × E isa EmptySet && E × B isa EmptySet
# test ∅ alias
@test B × E isa EmptySet

cpa = CartesianProductArray([B, 2.*B, 3.*B])
@test cpa * E isa EmptySet && E * cpa isa EmptySet
@test cpa × E isa EmptySet && E × cpa isa EmptySet

# testing cp of empty set with itself
@test ∅ * ∅ == ∅

# testing that the empty set is an absorbing element for the Minkowski sum
@test B + E isa EmptySet && E + B isa EmptySet
# testing the mathematical alias ⊕ 
@test B ⊕ E isa EmptySet && E ⊕ B isa EmptySet

msa = MinkowskiSumArray([B, 2.*B, 3.*B])
@test msa + E isa EmptySet && E + msa isa EmptySet
@test msa ⊕ E isa EmptySet && E ⊕ msa isa EmptySet

# testing M-sum of empty set with itself
@test ∅ + ∅ == ∅

# testing that the emptyset is neutral for the convex hull
@test CH(B, ∅) == B
@test CH(∅, B) == B

# test convex hull of empty set with itself
@test CH(∅, ∅) == ∅

# an_element function
@test_throws ErrorException an_element(∅)
