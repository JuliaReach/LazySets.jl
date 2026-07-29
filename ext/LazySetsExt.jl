# NOTE: This file and all included sub-files are organized like a package
#       extension for LazySets. While Julia itself is happy with that, the
#       package registrator does not accept a recursive package definition.
#       Hence, we simply include these files in the main LazySets code.

module LazySetsExt

import LazySets

include("LazySets/Ball1Ext.jl")
include("LazySets/BallpExt.jl")
include("LazySets/EmptySetExt.jl")
include("LazySets/HalfSpaceExt.jl")
include("LazySets/HParallelotopeExt.jl")
include("LazySets/HPolyhedronExt.jl")
include("LazySets/HPolytopeExt.jl")
include("LazySets/HyperplaneExt.jl")
include("LazySets/IntervalExt.jl")
include("LazySets/LineExt.jl")
include("LazySets/Line2DExt.jl")
include("LazySets/LineSegmentExt.jl")
include("LazySets/PolygonExt.jl")
include("LazySets/SparsePolynomialZonotopeExt.jl")
include("LazySets/StarExt.jl")
include("LazySets/TetrahedronExt.jl")
include("LazySets/UniverseExt.jl")
include("LazySets/VPolygonExt.jl")
include("LazySets/VPolytopeExt.jl")
include("LazySets/ZeroSetExt.jl")

end  # module
