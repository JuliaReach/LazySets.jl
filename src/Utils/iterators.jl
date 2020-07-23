# ========================
# iterator with no element
# ========================

struct EmptyIterator{T}
end

function Base.length(::EmptyIterator)
    return 0
end

function Base.eltype(::Type{EmptyIterator{T}}) where {T}
    return T
end

function Base.iterate(::EmptyIterator, state=nothing)
    return nothing
end

# ======================
# iterator over a vector
# ======================

struct VectorIterator{T}
    vector::T
end

function Base.length(it::VectorIterator)
    return length(it.vector)
end


function Base.eltype(::Type{VectorIterator{T}}) where {T}
    return eltype(T)
end

function Base.iterate(it::VectorIterator, state=1)
    if state > length(it)
        return nothing
    end
    element = it.vector[state]
    state += 1
    return (element, state)
end
