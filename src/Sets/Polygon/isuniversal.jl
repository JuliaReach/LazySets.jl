function isuniversal(P::Polygon, witness::Bool=false)
    # TODO support witness generation
    # return witness ? (false, non_element(P)) : false
    if witness
        throw(ArgumentError("witness generation is currently not supported"))
    end
    return false
end
