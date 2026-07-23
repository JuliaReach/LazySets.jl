module ExponentialUtilitiesExt

import ExponentialUtilities
using ExponentialUtilities: expv
using LazySets: exponential_backend, set_exponential_backend!
import LazySets: _expmv

function __init__()
    if ismissing(exponential_backend[])
        set_exponential_backend!(ExponentialUtilities)
    end
end

function _expmv(::Val{:ExponentialUtilities}, t, A, b; kwargs...)
    return expv(t, A, b; kwargs...)
end

end  # module
