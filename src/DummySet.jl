export DummySet

import Base:+,*

"""
    DummySet <: LazySet

Type that represents a "neutral" set with special properties, in a sense described
below.

### Fields

- `dim` -- ambient dimension of the `DummySet`
"""
struct DummySet <: LazySet
    dim::Int64
end

#=
# DummySet is the neutral element for +
# NOTE: We do not check for dimension mismatch.

function +(::DummySet, x)
    return x
end

function +(x, ::DummySet)
    return x
end

function +(::DummySet, x::LazySet)
    return x
end

function +(x::LazySet, ::DummySet)
    return x
end

function +(v::DummySet, ::DummySet)
    return v
end

# DummySet is the absorbing element for *
# NOTE: We do not check for dimension mismatch.
function *(x, v::DummySet)
    return v
end

function *(v::DummySet, x)
    return v
end

function *(::LazySet, v::DummySet)
    return v
end

function *(v::DummySet, ::LazySet)
    return v
end

function *(v::DummySet, ::DummySet)
    # NOTE: The result has the dimension of the first argument.
    # This is an arbitrary decision in the interest of avoiding additional
    # memory waste.
    return v
end

# dimension of a DummySet
function dim(V::DummySet)::Int64
    return V.dim
end

# support vector of a DummySet
function σ(d::AbstractVector{Float64}, V::DummySet)::Vector{Float64}
    #error("evaluating the support vector of a void set is undefined")
    return zeros(length(d))
end
=#
