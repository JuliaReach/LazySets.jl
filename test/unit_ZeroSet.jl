Z = ZeroSet(2)
B = BallInf(ones(2), 1.)

# testing that the zero set is neutral element for the Minkowski sum
@test B ⊕ Z == B && Z ⊕ B == B

cpa = MinkowskiSumArray([B, 2.*B, 3.*B])
@test cpa ⊕ Z == cpa && Z ⊕ cpa == cpa

# test M-sum of zero set with itself
@test ZeroSet(2) ⊕ ZeroSet(2) == ZeroSet(2)

# an_element function
@test an_element(Z) ∈ Z

# subset
z = ZeroSet(1)
s1 = Singleton([0.])
s2 = Singleton([2.])
@test ⊆(z, s1) && ⊆(z, s1, true)[1]
subset, point = ⊆(z, s2, true)
@test !⊆(z, s2) && !subset && point ∈ z && !(point ∈ s2)
@test ⊆(z, z) && !⊆(z, ZeroSet(2))
