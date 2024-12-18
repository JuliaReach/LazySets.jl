function isuniversal(U::Universe{N}, witness::Bool=false) where {N}
    return witness ? (true, N[]) : true
end
