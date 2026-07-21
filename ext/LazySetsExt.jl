# NOTE: This file and all included sub-files are organized like a package
#       extension for LazySets. While Julia itself is happy with that, the
#       package registrator does not accept a recursive package definition.
#       Hence, we simply include these files in the main LazySets code.

module LazySetsExt

import LazySets

include("LazySets/LazySetsBall1Ext.jl")
include("LazySets/LazySetsBallpExt.jl")
include("LazySets/LazySetsEmptySetExt.jl")
include("LazySets/LazySetsHalfSpaceExt.jl")
include("LazySets/LazySetsHParallelotopeExt.jl")
include("LazySets/LazySetsHPolyhedronExt.jl")
include("LazySets/LazySetsHPolytopeExt.jl")
include("LazySets/LazySetsHyperplaneExt.jl")
include("LazySets/LazySetsIntervalExt.jl")
include("LazySets/LazySetsLineExt.jl")
include("LazySets/LazySetsLine2DExt.jl")
include("LazySets/LazySetsLineSegmentExt.jl")
include("LazySets/LazySetsPolygonExt.jl")
include("LazySets/LazySetsSparsePolynomialZonotopeExt.jl")
include("LazySets/LazySetsStarExt.jl")
include("LazySets/LazySetsTetrahedronExt.jl")
include("LazySets/LazySetsUniverseExt.jl")
include("LazySets/LazySetsVPolygonExt.jl")
include("LazySets/LazySetsVPolytopeExt.jl")
include("LazySets/LazySetsZeroSetExt.jl")

end  # module
