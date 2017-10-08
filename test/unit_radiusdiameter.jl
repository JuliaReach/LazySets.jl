# Radius and Diameter Approximation of a centered unit Ball2
b = Ball2([0., 0., 0.], 1.)
rexp = norm([1.,1.,1.])
dexp = norm([2.,2.,2.])
@test radius_approximation(b) >= rexp
@test diameter_approximation(b) >= dexp

# Radius and Diameter Approximation of a not-centered unit Ball2
b = Ball2([1., 2., 0.], 1.)
rexp = norm([2.,3.,1.])
dexp = norm([2.,2.,2.])
@test radius_approximation(b) >= rexp
@test diameter_approximation(b) >= dexp

# Radius and Diameter Approximation of a non-unit Ball2
b = Ball2([0., 0., 0.], 3.)
rexp = norm([3.,3.,3.])
dexp = norm([6.,6.,6.])
@test radius_approximation(b) >= rexp
@test diameter_approximation(b) >= dexp