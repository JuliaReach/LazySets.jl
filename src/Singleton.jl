export Singleton

"""
    Singleton{N<:Real} <: AbstractSingleton{N}

Type that represents a singleton, that is, a set with a unique element.

### Fields

- `element` -- the only element of the set
"""
struct Singleton{N<:Real} <: AbstractSingleton{N}
    element::Vector{N}
end


# --- AbstractSingleton interface functions ---


"""
    element(S::Singleton{N})::Vector{N} where {N<:Real}

Return the element of a singleton.

### Input

- `S` -- singleton

### Output

The element of the singleton.
"""
function element(S::Singleton{N})::Vector{N} where {N<:Real}
    return S.element
end

"""
    element(S::Singleton{N}, i::Int)::N where {N<:Real}

Return the i-th entry of the element of a singleton.

### Input

- `S` -- singleton
- `i` -- dimension

### Output

The i-th entry of the element of the singleton.
"""
function element(S::Singleton{N}, i::Int)::N where {N<:Real}
    return S.element[i]
end
