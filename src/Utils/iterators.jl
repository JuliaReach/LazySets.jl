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
