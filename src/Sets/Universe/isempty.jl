function isempty(U::Universe, witness::Bool=false)
    if witness
        return (false, an_element(U))
    end
    return false
end
