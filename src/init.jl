function __init__()
    @require CDDLib = "3391f64e-dcde-5f30-b752-e11513730f60" include("Initialization/init_CDDLib.jl")
    @require ExponentialUtilities = "d4d017d3-3776-5f7e-afef-a10c40355c18" include("Initialization/init_ExponentialUtilities.jl")
    @require IntervalMatrices = "5c1f47dc-42dd-5697-8aaa-4d102d140ba9" include("Initialization/init_IntervalMatrices.jl")
    @require IntervalBoxes = "43d83c95-ebbb-40ec-8188-24586a1458ed" include("Initialization/init_IntervalBoxes.jl")
    @require MiniQhull = "978d7f02-9e05-4691-894f-ae31a51d76ca" include("Initialization/init_MiniQhull.jl")
    @require Optim = "429524aa-4258-5aef-a3af-852621145aeb" include("Initialization/init_Optim.jl")
    @require Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029" begin
        include("Initialization/init_Polyhedra.jl")
        @require GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326" include("Initialization/init_GeometryBasics.jl")
    end
    @require RangeEnclosures = "1b4d18b6-9e5d-11e9-236c-f792b01831f8" include("Initialization/init_RangeEnclosures.jl")
    @require StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c" include("Initialization/init_StaticArraysCore.jl")
    @require TaylorModels = "314ce334-5f6e-57ae-acf6-00b6e903104a" include("Initialization/init_TaylorModels.jl")

    return nothing
end
