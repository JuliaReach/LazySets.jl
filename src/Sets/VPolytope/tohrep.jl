# see ext/LazySets/LazySetsVPolytopeExt.jl

# see ext/Polyhedra/PolyhedraVPolytopeExt.jl
function _tohrep(P; backend)
    mod = Base.get_extension(@__MODULE__, :LazySetsPolyhedraExt)
    require(mod, :Polyhedra; fun_name="tohrep")
    error()
end
