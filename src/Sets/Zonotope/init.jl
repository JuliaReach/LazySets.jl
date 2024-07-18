function __init__()
    @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" include("init_LazySets.jl")
    @require StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c" include("init_StaticArraysCore.jl")
end
