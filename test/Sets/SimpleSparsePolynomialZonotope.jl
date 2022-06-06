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

    @test overapproximate(S, Zonotope) == Z
end
