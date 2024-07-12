function polynomial_order(P::SSPZ)
    return maximum(sum, eachcol(expmat(P)))
end
