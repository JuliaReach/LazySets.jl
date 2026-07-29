# currently supported packages
const SUPPORTED_SDP_PACKAGES = [:SCS]

# =============
# Global option
# =============

const sdp_solver = Ref{Any}(missing)  # global state of the SDP solver

function set_sdp_solver!(solver::Module)
    sdp_solver[] = Val(Symbol(solver))
    return sdp_solver[]
end

function get_sdp_solver()
    if ismissing(sdp_solver[])
        throw(ArgumentError("no semidefinite-programming (SDP) solver is loaded; load one of " *
                            "these packages: $SUPPORTED_SDP_PACKAGES"))
    end
    return sdp_solver[]
end

# default semidefinite-programming solver
function default_sdp_solver()
    solver = get_sdp_solver()
    return _default_sdp_solver(solver)
end

# see ext/LazySetsSCSExt.jl
function _default_sdp_solver(solver)
    error()
end
