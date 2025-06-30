for N in [Float64, Float32, Rational{Int}]
    # an_element default implementation
    U = Universe{N}(2)
    @test_throws ArgumentError LazySets._an_element_lazySet(U)
end
