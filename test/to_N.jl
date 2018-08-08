# numbers

function to_N(N::Type{NUM}, v::NUM) where {NUM<:Real}
    return v
end

function to_N(N::Type{<:Real}, v::Real)
    return N(v)
end

function to_N(N::Type{<:Rational{NUM}}, v::Real) where {NUM}
    return rationalize(NUM, v)
end

# arrays

function to_N(N::Type{NUM}, a::AbstractVector{NUM}) where {NUM<:Real}
    return a
end

function to_N(N::Type{<:Real}, a::AbstractVector{<:Real})
    r = Vector{N}(undef, length(a))
    for i in eachindex(a)
        r[i] = to_N(N, a[i])
    end
    return r
end

# matrices

function to_N(N::Type{NUM}, m::AbstractMatrix{NUM}) where {NUM<:Real}
    return m
end

function to_N(N::Type{<:Real}, m::AbstractMatrix{<:Real})
    r = Matrix{N}(undef, size(m))
    for j = 1:size(m,2)
        for i = 1:size(m,1)
            r[i, j] = to_N(N, m[i, j])
        end
    end
    return r
end

# arrays of arrays

function to_N(N::Type{NUM}, a::AbstractVector{<:AbstractVector{NUM}}) where {NUM<:Real}
    return a
end

function to_N(N::Type{<:Real}, a::AbstractVector{<:AbstractVector{<:Real}})
    r = Vector{Vector{N}}(undef, length(a))
    for i in eachindex(a)
        v = a[i]
        b = Vector{N}(undef, length(v))
        for j in eachindex(v)
            b[j] = to_N(N, v[j])
        end
        r[i] = b
    end
    return r
end
