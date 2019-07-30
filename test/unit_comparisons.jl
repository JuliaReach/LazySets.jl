using LazySets: _leq, _geq, isapproxzero, _isapprox, _ztol

# approximate <= and
@test _leq(2e-15, 1e-15) && _leq(1e-15, 2e-15)
@test _leq(1//100, 1//99) && _leq(1//100, 0.099)

# approximate <= with ztol
@test !_leq(2e-15, 1e-15, ztol=1e-15) && _leq(1e-15, 2e-15, ztol=1e-15)

# approximate >=
@test _geq(2e-15, 1e-15) && _geq(1e-15, 2e-15)

# approximate >= with ztol
@test _geq(2e-15, 1e-15, ztol=1e-15) && !_geq(1e-15, 2e-15, ztol=1e-15)

# default absolute zero tolerance for rational
_ztol(eltype(1/100)) == zero(Rational{Int})

# default absolute zero tolerance for FP
_ztol(eltype(0.01)) == sqrt(eps(Float64))

# approximately zero tests
@test isapproxzero(0//1)
@test isapproxzero(1e-8) && !isapproxzero(1e-8, ztol=1e-10)

# approximate equality
@test !_isapprox(1//1, 1//1 + 1//1000000000000)
@test _isapprox(2e-15, 1e-15)
@test !_isapprox(2e-6, 1e-6)
