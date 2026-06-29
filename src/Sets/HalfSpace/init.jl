function __init__()
    @require SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8" include("init_SymEngine.jl")
    @require Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7" include("init_Symbolics.jl")
end
