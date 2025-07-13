using LazySets, BenchmarkTools, Profile

MZ = rand(MatrixZonotope; dim=(5, 5), num_generators=10)
MZ2 = rand(MatrixZonotope; dim=(5, 5), num_generators=10)
S = rand(SparsePolynomialZonotope; dim=5, nparams=3, maxdeg=3,
         num_dependent_generators=8,
         num_independent_generators=5)

lm = MZ * S
#@btime overapproximate(MZ * S)
exp = ExponentialMap(MZ*MZ2, S)
@profile overapproximate(exp, 2)