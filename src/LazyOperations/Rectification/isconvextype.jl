function isconvextype(::Type{<:Rectification})
    return false
end

function isconvextype(::Type{<:Rectification{N,<:AbstractHyperrectangle}}) where {N}
    return true
end

function isconvextype(::Type{<:Rectification{N,<:EmptySet}}) where {N}
    return true
end

function isconvextype(::Type{<:Rectification{N,<:Universe}}) where {N}
    return true
end
