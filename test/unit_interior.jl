P = BallInf(zeros(2), 0.15)
d = [0.15,0.15]
@test !interior(d, P)
d = [0.15,0.15] .- LazySets._TOL_F64.rtol
@test interior(d, P)
