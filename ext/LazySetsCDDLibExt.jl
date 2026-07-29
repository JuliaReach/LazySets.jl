module LazySetsCDDLibExt

using CDDLib: Library  # NOTE: this is an internal function
import LazySets: default_cddlib_backend

function default_cddlib_backend(::Type{<:AbstractFloat})
    return Library()
end

function default_cddlib_backend(::Type{<:Rational})
    return Library(:exact)
end

end  # module
