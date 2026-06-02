@validate function translate(msa::MinkowskiSumArray, x::AbstractVector)
    return MinkowskiSumArray([translate(X, x) for X in array(msa)])
end
