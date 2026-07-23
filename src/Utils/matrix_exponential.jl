# currently supported packages
const SUPPORTED_EXPONENTIAL_PACKAGES = [:ExponentialUtilities, :Expokit]

# =============
# Global option
# =============

const exponential_backend = Ref{Any}(missing)  # global state of the exponential backend

function set_exponential_backend!(backend::Module)
    exponential_backend[] = Val(Symbol(backend))
    return exponential_backend[]
end

function get_exponential_backend()
    if ismissing(exponential_backend[])
        throw(ArgumentError("no exponential backend is loaded; load one of these packages: " *
                            "$SUPPORTED_EXPONENTIAL_PACKAGES"))
    end
    return exponential_backend[]
end

# use default backend
function _expmv(t, A, b; kwargs...)
    backend = get_exponential_backend()
    return _expmv(backend, t, A, b; kwargs...)
end

# see ext/ExpokitExt.jl and ext/ExponentialUtilitiesExt.jl
function _expmv(backend, t, A, b; kwargs...)
    # note: `@__MODULE__` does not detect when a package extension is loaded;
    #       however, the method will simply be overwritten
    require(@__MODULE__, SUPPORTED_EXPONENTIAL_PACKAGES; fun_name="_expmv")
    error()
end
