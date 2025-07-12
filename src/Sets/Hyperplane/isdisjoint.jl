@validate function isdisjoint(hp1::Hyperplane, hp2::Hyperplane, witness::Bool=false)
    return _isdisjoint_hyperplane_hyperplane(hp1, hp2, witness)
end

function _isdisjoint_hyperplane_hyperplane(hp1, hp2, witness::Bool=false)
    require(@__MODULE__, :LazySets; fun_name="isdisjoint")

    if isequivalent(hp1, hp2)
        res = false
        if witness
            w = an_element(hp1)
        end
    else
        cap = intersection(hp1, hp2)
        res = cap isa EmptySet
        if !res && witness
            w = an_element(cap)
        end
    end
    if res
        return _witness_result_empty(witness, true, hp1, hp2)
    end
    return witness ? (false, w) : false
end
