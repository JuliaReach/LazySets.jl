#=
This file defines internal functions for compatibility across
different Julia versions.
=#

using Compat
using Compat: copyto!, axes, argmax, @warn
import Compat.String
using Compat.LinearAlgebra
import Compat.LinearAlgebra: norm, checksquare, LAPACKException,
                             SingularException, ×, cond
import Compat.InteractiveUtils.subtypes
export _At_mul_B

@static if VERSION < v"0.7-"
    using Compat.Random
    import Compat.Random: GLOBAL_RNG
    @inline _At_mul_B(A, B) = At_mul_B(A, B)
    @inline _At_ldiv_B(A, B) = At_ldiv_B(A, B)
    expmat = expm
else
    using SparseArrays
    using Random
    using Random: GLOBAL_RNG, SamplerType
    @inline _At_mul_B(A, B) = transpose(A) * B
    @inline _At_ldiv_B(A, B) = transpose(A) \ B
    expmat = exp
end
