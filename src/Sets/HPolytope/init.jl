function __init__()
    @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" include("init_LazySets.jl")
    @require Polyhedra = "67491407-f73d-577b-9b50-8179a7c68029" eval(load_polyhedra_hpolytope())
    @require Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7" eval(load_symbolics_hpolytope())
end
