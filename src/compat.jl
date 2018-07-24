#=
This file defines internal functions for compatibility across
different Julia versions.
=#

using Compat
import Compat.String

if VERSION < v"0.7-"
    import Base.LinAlg:norm, checksquare
    import Base: eye, ×
else
    using SparseArrays, LinearAlgebra
    import LinearAlgebra:norm, checksquare, eye, ×
end

export _At_mul_B

@inline function _At_mul_B(A, B)
    if VERSION > v"0.7-"
        transpose(A) * B
    else
        At_mul_B(A, B)
    end
end
