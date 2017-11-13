import LazySets.Approximations:norm,
                               radius,
                               diameter

# =======
#  Ball2
# =======

# approximation of a centered unit Ball2
b = Ball2([0., 0., 0.], 1.)
@test norm(b) ≈ 1.  # in the infinity norm (default)
@test radius(b) ≈ 1.  # in the infinity norm (default)
@test diameter(b) ≈ 2.  # in the infinity norm (default)

# approximations of a non-centered unit Ball2
b = Ball2([1., 2., 0.], 1.)
@test norm(b) ≈ 3.  # in the infinity norm (default)
@test radius(b) ≈ 1.  # in the infinity norm (default)
@test diameter(b) ≈ 2.  # in the infinity norm (default)


# ================
#  Hyperrectangle
# ================

# metrics in the infinity norm (default)
h = Hyperrectangle([-1., 2.], [0.2, 0.5])
@test norm(h, Inf) ≈ 2.5
@test radius(h, Inf) ≈ 0.5
@test diameter(h, Inf) ≈ 1.

# metrics in the 2-norm
@test norm(h, 2) ≈ sqrt(2.5^2 + 1.2^2)
@test radius(h, 2) ≈ sqrt(.5^2 + .2^2)
@test diameter(h, 2) ≈ 2*sqrt(.5^2 + .2^2)

# =========
#  BallInf
# =========

# metrics in the infinity norm (default)
b = BallInf([-1., -2., -4], 0.2)
@test norm(b, Inf) ≈ 4.2
@test radius(b, Inf) ≈ 0.2
@test diameter(b, Inf) ≈ 0.4

# metrics in the 2-norm
@test norm(b, 2) ≈ sqrt(1.2^2 + 2.2^2 + 4.2^2)
@test radius(b, 2) ≈ sqrt(0.2^2 * 3)
@test diameter(b, 2) ≈ 2*sqrt(0.2^2 * 3)
