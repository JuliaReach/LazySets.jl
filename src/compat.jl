#=
This file defines internal functions for compatibility across
different Julia versions.
=#

export _At_mul_B

function _At_mul_B(A, B)
    if VERSION > v"0.7-"
        transpose(A) * B
    else
        At_mul_B(A, B)
    end
end
