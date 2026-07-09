using LazySets, Test

for N in [Float64, Float32]
    # Test Kamenev on a 2D unit ball in the 2-norm
    B = Ball2(zeros(N, 2), one(N))
    tol = N(1e-2)
    inner, outer, err = kamenev(B, tol)
    @test inner ⊂ B
    @test B ⊂ outer
   
    # Test underapproximate/overapproximate aliases
    inner2 = underapproximate(B, VPolytope, tol)
    @test inner2 == inner

    outer2 = overapproximate(B, HPolytope, tol)
    @test outer2 == outer
end
