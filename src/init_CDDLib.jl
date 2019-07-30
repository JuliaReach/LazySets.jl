default_cddlib_backend(::Type{<:AbstractFloat}) = CDDLib.Library()
default_cddlib_backend(::Type{<:Rational}) = CDDLib.Library(:Exact)
