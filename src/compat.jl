#=
This file defines internal functions for compatibility across
different Julia versions.
=#

using Compat
using Compat: copyto!, axes
import Compat.String
export _At_mul_B

if VERSION < v"0.7-"
    import Base.LinAlg: norm, checksquare, LAPACKException, SingularException
    import Base: eye, ×
    @inline function _At_mul_B(A, B)
        return At_mul_B(A, B)
    end
else
    using SparseArrays, LinearAlgebra
    import LinearAlgebra: norm, checksquare, LAPACKException, SingularException,
                          eye, ×
    @inline function _At_mul_B(A, B)
        return transpose(A) * B
    end
end
