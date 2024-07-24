# _is_linearcombination
@test LazySets._is_linearcombination(:(2*x1 - 4))
@test LazySets._is_linearcombination(:(6.1 - 5.3*f - 0.1*g))
@test !LazySets._is_linearcombination(:(2*x1^2))
@test !LazySets._is_linearcombination(:(x1^2 - 4*x2 + x3 + 2))

# _free_symbols
@test LazySets._free_symbols(:(x1 = 1)) == [SymEngine.Basic(:(x1))]
@test LazySets._free_symbols(:(2*x2 <= 4)) == [SymEngine.Basic(:(x2))]
@test_throws ErrorException LazySets._free_symbols(:(x3 != 4))
