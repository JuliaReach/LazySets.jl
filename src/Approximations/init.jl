function __init__()
    @require ExponentialUtilities = "d4d017d3-3776-5f7e-afef-a10c40355c18" (nothing,)
    @require IntervalMatrices = "5c1f47dc-42dd-5697-8aaa-4d102d140ba9" include("init_IntervalMatrices.jl")
    @require StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c" include("init_StaticArraysCore.jl")
    @require TaylorModels = "314ce334-5f6e-57ae-acf6-00b6e903104a" include("init_TaylorModels.jl")
end
