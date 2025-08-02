using LazySets, Test

@static if isdefined(@__MODULE__, :CDDLib)
    @test LazySets.default_cddlib_backend(Float64) == CDDLib.Library(:float)
    @test LazySets.default_cddlib_backend(Float32) == CDDLib.Library(:float)
    @test LazySets.default_cddlib_backend(Rational{Int}) == CDDLib.Library(:exact)
end
