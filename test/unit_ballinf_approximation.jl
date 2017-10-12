# BallInf approximation of a 3D unit ball centered at [1,2,0] of radius 1
b = Ball2([1., 2., 0.], 1.)
bi = ballinf_approximation(b)
biexp = BallInf([1., 2., 0.], 1.)
@test bi.center ≈ biexp.center
@test bi.radius ≈ biexp.radius

# BallInf approximation of a 2D polygon whose vertices are
# (0,1), (1,1), (2,-1), (-1,0)
p = HPolygon()
addconstraint!(p, LinearConstraint([0., 1.], 1.))
addconstraint!(p, LinearConstraint([-1., 1.], 1.))
addconstraint!(p, LinearConstraint([-1., -3.], 1.))
addconstraint!(p, LinearConstraint([2., 1.], 3.))
bi = ballinf_approximation(p)
biexp = BallInf([.5, 0.], 1.5)
@test bi.center ≈ biexp.center
@test bi.radius ≈ biexp.radius
