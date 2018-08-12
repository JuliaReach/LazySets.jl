#=
This file defines internal functions for compatibility across
different Julia versions.
=#

using Compat
using Compat: copyto!, axes, argmax
import Compat.String
using Compat.LinearAlgebra
import Compat.LinearAlgebra: norm, checksquare, LAPACKException,
                             SingularException, eye, ×
import Compat.InteractiveUtils.subtypes
export _At_mul_B

if VERSION < v"0.7-"
    @inline function _At_mul_B(A, B)
        return At_mul_B(A, B)
    end
else
    using SparseArrays
    @inline function _At_mul_B(A, B)
        return transpose(A) * B
    end
end
