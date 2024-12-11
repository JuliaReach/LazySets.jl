function isuniversal(∅::EmptySet{N}, witness::Bool=false) where {N}
    if witness
        return (false, zeros(N, dim(∅)))
    else
        return false
    end
end
