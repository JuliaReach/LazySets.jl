P = BallInf(zeros(2), 0.15)
d = [0.15,0.15]
@test !is_interior_point(d, P)
@test !is_interior_point(d, P; p=2.0)
@test !is_interior_point(d, P; ε=1.0)
@test is_interior_point(d, P; ε=0.0)

d = [0.15,0.15] .- LazySets._TOL_F64.rtol
@test is_interior_point(d, P)
@test is_interior_point(d, P; p=2.0)
