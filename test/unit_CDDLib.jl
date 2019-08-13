using LazySets: default_cddlib_backend

@test default_cddlib_backend(Float64) == CDDLib.Library(:float)
@test default_cddlib_backend(Float32) == CDDLib.Library(:float)
@test default_cddlib_backend(Rational{Int}) == CDDLib.Library(:exact)
