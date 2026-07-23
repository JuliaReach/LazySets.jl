function __init__()
    @require IntervalMatrices = "5c1f47dc-42dd-5697-8aaa-4d102d140ba9" include("Initialization/init_IntervalMatrices.jl")
    @require Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029" begin
        include("Initialization/init_Polyhedra.jl")
        @require GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326" include("Initialization/init_GeometryBasics.jl")
    end
    @require StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c" include("Initialization/init_StaticArraysCore.jl")
    @require TaylorModels = "314ce334-5f6e-57ae-acf6-00b6e903104a" include("Initialization/init_TaylorModels.jl")

    return nothing
end
