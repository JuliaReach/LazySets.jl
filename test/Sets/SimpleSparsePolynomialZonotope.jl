@testset "Simple Sparse Polynomial Zonotope" begin

    # example from slide 13 of Niklas talk at JRJID3
    c = [2.0, 0.0]
    G = [1 2;2 2.]
    E = [1 4;1 2]

    S = SimpleSparsePolynomialZonotope(c, G, E)

    @test genmat(S) == G
    @test dim(S) == 2
    @test ngens(S) == 2
    @test nparams(S) == 2
    @test order(S) == 1 // 1

    @test overapproximate(S, Zonotope) == Zonotope([3., 1], [1 1;2 1.])
    @test length(overapproximate(S, Zonotope; nsdiv=3)) == 9
    @test length(overapproximate(S, Zonotope; partition=(2, 3))) == 6
end
