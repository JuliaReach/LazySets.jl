module LazySetsSCSExt

using LazySets: sdp_solver, set_sdp_solver!
import SCS
using SCS: Optimizer
using MathOptInterface: OptimizerWithAttributes, Silent
import LazySets: _default_sdp_solver

function __init__()
    if ismissing(sdp_solver[])
        set_sdp_solver!(SCS)
    end
end

function _default_sdp_solver(::Val{:SCS})
    return OptimizerWithAttributes(Optimizer, Silent() => true)
end

end  # module
