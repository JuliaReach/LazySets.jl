function __init__()
    @require Expokit = "a1e7a1ef-7a5d-5822-a38c-be74e1bb89f4" include("init_Expokit.jl")
    @require IntervalMatrices = "5c1f47dc-42dd-5697-8aaa-4d102d140ba9" include("init_IntervalMatrices.jl")
    @require IntervalConstraintProgramming = "138f1668-1576-5ad7-91b9-7425abbf3153" include("init_IntervalConstraintProgramming.jl")
    @require Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9" include("init_Ipopt.jl")
    @require StaticArrays = "90137ffa-7385-5640-81b9-e52037218182" include("init_StaticArrays.jl")
    @require TaylorModels = "314ce334-5f6e-57ae-acf6-00b6e903104a" include("init_TaylorModels.jl")
end
