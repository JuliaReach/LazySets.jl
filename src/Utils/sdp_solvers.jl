# currently supported packages
const SUPPORTED_SDP_PACKAGES = [:SCS]

# =============
# Global option
# =============

global sdp_solver = missing  # global state of the SDP solver

function set_sdp_solver!(solver::Module)
    return global sdp_solver = Val(Symbol(solver))
end

function get_sdp_solver()
    ismissing(sdp_solver) && error("no semi-definite-SDP solver is loaded; load one of these " *
                                   "packages: $SUPPORTED_SDP_PACKAGES")
    return sdp_solver
end

# default semidefinite-programming solver
function default_sdp_solver()
    solver = get_sdp_solver()
    return _default_sdp_solver(solver)
end

# ===
# SCS
# ===

function load_scs()
    return quote
        if ismissing(sdp_solver)
            set_sdp_solver!(SCS)
        end
    end
end  # quote / load_scs

function _default_sdp_solver(::Val{:SCS})
    solver = JuMP.optimizer_with_attributes(() -> SCS.Optimizer())
    JuMP.set_attribute(solver, JuMP.MOI.Silent(), true)
    return solver
end
