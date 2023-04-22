# currently supported packages
const SUPPORTED_EXPONENTIAL_PACKAGES = [:ExponentialUtilities, :Expokit]

# =============
# Global option
# =============

global exponential_backend = missing  # global state of the exponential backend

function set_exponential_backend!(backend::Module)
    global exponential_backend = Val(Symbol(backend))
end

function get_exponential_backend()
    ismissing(exponential_backend) && error("no exponential backend is loaded; " *
        "load one of these packages: $SUPPORTED_EXPONENTIAL_PACKAGES")
    return exponential_backend
end

# use default backend
function _expmv(t, A, b; kwargs...)
    backend = get_exponential_backend()
    return _expmv(backend, t, A, b; kwargs...)
end

# ====================
# ExponentialUtilities
# ====================

function load_exponentialutilities()
return quote

if ismissing(exponential_backend)
    set_exponential_backend!(ExponentialUtilities)
end

end end  # quote / load_exponentialutilities

function _expmv(::Val{:ExponentialUtilities}, t, A, b; kwargs...)
    return ExponentialUtilities.expv(t, A, b; kwargs...)
end

# =======
# Expokit
# =======

function load_expokit()
return quote

if ismissing(exponential_backend)
    set_exponential_backend!(Expokit)
end

end end  # quote / load_expokit

function _expmv(::Val{:Expokit}, t, A, b::Vector; kwargs...)
    return Expokit.expmv(t, A, b; kwargs...)
end

function _expmv(backend::Val{:Expokit}, t, A, b; kwargs...)
    b = Vector(b)  # Expokit requires a dense vector
    return _expmv(backend, t, A, b; kwargs...)
end
