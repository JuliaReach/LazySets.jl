using LazySets, Test
using LazySets.ReachabilityBase.Basetype
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    Z = ZeroSet{N}(2)

    # basetype
    @test basetype(Z) == basetype(typeof(Z)) == ZeroSet

    # eltype
    @test eltype(Z) == eltype(typeof(Z)) == N

    # an_element default implementation
    U = Universe{N}(2)
    @test_throws ArgumentError LazySets._an_element_lazyset(U)
end
