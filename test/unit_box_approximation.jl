# ==============================
# Testing box approximation
# ==============================
# Box approximation of a 2D square
b = BallInf([1., 1.], 0.1)
h = box_approximation(b)
hexp = Hyperrectangle([1., 1.], [0.1, 0.1])
@test h.center ≈ hexp.center
@test h.radius ≈ hexp.radius

# Box approximation of a 2D unit ball
b = Ball2([1., -2.], 0.2)
h = box_approximation(b)
hexp = Hyperrectangle([1., -2.], [0.2, 0.2])
@test h.center ≈ hexp.center
@test h.radius ≈ hexp.radius


# Box approximation of a 3D unit ball
b = Ball2([1., 2., 0.], 1.)
h = box_approximation(b)
hexp = Hyperrectangle([1., 2., 0.], [1., 1., 1.])
@test h.center ≈ hexp.center
@test h.radius ≈ hexp.radius


# ===================================================================
# Testing box_approximation_symmetric (= symmetric interval hull)
# ===================================================================
# Box approximation of a 2D square
b = BallInf([1., 1.], 0.1)
h = box_approximation_symmetric(b)
hexp = Hyperrectangle([0.0, 0.0], [1.1, 1.1])
@test h.center ≈ hexp.center
@test h.radius ≈ hexp.radius

# Box approximation of a 2D unit ball
b = Ball2([1., -2.], 0.2)
h = box_approximation_symmetric(b)
hexp = Hyperrectangle([0., 0.], [1.2, 2.2])
@test h.center ≈ hexp.center
@test h.radius ≈ hexp.radius

# Box approximation of a 3D unit ball
b = Ball2([1., 2., 0.], 0.1)
h = box_approximation_symmetric(b)
hexp = Hyperrectangle([0., 0., 0.], [1.1, 2.1, 0.1])
@test h.center ≈ hexp.center
@test h.radius ≈ hexp.radius

# Box approximation of a 4D hyperrectangle
b = Hyperrectangle([-1.5, -2.5, 2.4, -0.4], [0.1, 0.2, 0.3, 0.4])
h = box_approximation_symmetric(b)
hexp = Hyperrectangle(zeros(4), [1.6, 2.7, 2.7, 0.8])
@test h.center ≈ hexp.center
@test h.radius ≈ hexp.radius

# Testing alias symmetric_interval_hull
h = symmetric_interval_hull(b)

