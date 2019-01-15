import LazySets: _leq, _geq, isapproxzero, _isapprox

# approximate <= and >=
@test _leq(2e-15, 1e-15) && _leq(1e-15, 2e-15)
@test _geq(2e-15, 1e-15) && _geq(1e-15, 2e-15)

# approximate <= and >= with ztol
@test !_leq(2e-15, 1e-15, ztol=1e-15) && _leq(1e-15, 2e-15, ztol=1e-15)
@test _geq(2e-15, 1e-15, ztol=1e-15) && !_geq(1e-15, 2e-15, ztol=1e-15)

@test isapproxzero(0//1)
@test isapproxzero(1e-8) && !isapproxzero(1e-8, ztol=1e-10)

@test _isapprox(2e-15, 1e-15)
