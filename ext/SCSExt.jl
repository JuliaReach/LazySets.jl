module SCSExt

using LazySets: sdp_solver, set_sdp_solver!
import SCS
using SCS: MOI, Optimizer
import LazySets: _default_sdp_solver

function __init__()
    if ismissing(sdp_solver[])
        set_sdp_solver!(SCS)
    end
end

function _default_sdp_solver(::Val{:SCS})
    return MOI.OptimizerWithAttributes(Optimizer, MOI.Silent() => true)
end

end  # module
