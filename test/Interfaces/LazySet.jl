
using LazySets, Test
using LazySets.ReachabilityBase.Basetype

for N in [Float64, Float32, Rational{Int}]
    Z = ZeroSet{N}(2)

    # basetype
    @test basetype(Z) == basetype(typeof(Z)) == ZeroSet

    # eltype
    @test eltype(Z) == eltype(typeof(Z)) == N
end
