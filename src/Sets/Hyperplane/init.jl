function __init__()
    @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" include("init_LazySets.jl")
    @require SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8" include("init_SymEngine.jl")
    @require Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7" include("init_Symbolics.jl")
end
