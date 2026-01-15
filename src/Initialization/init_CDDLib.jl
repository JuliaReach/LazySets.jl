default_cddlib_backend(::Type{<:AbstractFloat}) = CDDLib.Library()  # NOTE: this is an internal function
default_cddlib_backend(::Type{<:Rational}) = CDDLib.Library(:exact)
