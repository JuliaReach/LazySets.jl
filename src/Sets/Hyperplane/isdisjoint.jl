@validate function isdisjoint(hp1::Hyperplane, hp2::Hyperplane, witness::Bool=false)
    return _isdisjoint_hyperplane_hyperplane(hp1, hp2, witness)
end

# see ext/LazySets/LazySetsHyperplaneExt.jl
_isdisjoint_hyperplane_hyperplane(hp1, hp2, witness=false) = error()
