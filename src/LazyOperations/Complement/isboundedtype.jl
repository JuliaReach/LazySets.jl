function isboundedtype(::Type{<:Complement})
    return false
end

function isboundedtype(::Type{<:Complement{N,<:Universe}}) where {N}
    return true
end
