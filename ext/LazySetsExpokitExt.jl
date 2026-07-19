module LazySetsExpokitExt

import Expokit
using Expokit: expmv
using LazySets: exponential_backend, set_exponential_backend!
import LazySets: _expmv

function __init__()
    if ismissing(exponential_backend[])
        set_exponential_backend!(Expokit)
    end
end

function _expmv(::Val{:Expokit}, t, A, b::Vector; kwargs...)
    return expmv(t, A, b; kwargs...)
end

function _expmv(backend::Val{:Expokit}, t, A, b; kwargs...)
    b = Vector(b)  # Expokit requires a dense vector
    return _expmv(backend, t, A, b; kwargs...)
end

end  # module
