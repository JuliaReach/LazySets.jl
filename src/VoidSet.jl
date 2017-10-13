export VoidSet

"""
    VoidSet <: LazySet

Type that represents a void (neutral) set with respect to Minkowski sum.

### Fields

- `dim` -- ambient dimension of the `VoidSet` 
"""
struct VoidSet <: LazySet
    dim::Int64
end

# VoidSet is the neutral element for +
# NOTE: We do not check for dimension mismatch.
import Base.+
function +(::VoidSet, x)
    return x
end

function +(x, ::VoidSet)
    return x
end

function +(::VoidSet, x::LazySet)
    return x
end

function +(x::LazySet, ::VoidSet)
    return x
end

function +(v::VoidSet, ::VoidSet)
    return v
end

# VoidSet is the absorbing element for *
# NOTE: We do not check for dimension mismatch.
import Base.*
function *(x, v::VoidSet)
    return v
end

function *(v::VoidSet, x)
    return v
end

function *(::LazySet, v::VoidSet)
    return v
end

function *(v::VoidSet, ::LazySet)
    return v
end

function *(v::VoidSet, ::VoidSet)
    # NOTE: The result has the dimension of the first argument.
    # This is an arbitrary decision in the interest of avoiding additional
    # memory waste.
    return v
end

# dimension of a VoidSet
function dim(V::VoidSet)::Int64
    return V.dim
end

# support vector of a VoidSet
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, V::VoidSet)::Vector{Float64}
    #error("evaluating the support vector of a void set is undefined")
    return zeros(length(d))
end
