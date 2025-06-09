function area(B::BallInf)
    n = dim(B)
    @assert n ∈ (2, 3) "this function only applies to two-dimensional or " *
                        "three-dimensional sets, but the given set is " *
                        "$n-dimensional"

    if n == 2
        return (2 * B.radius)^2
    else
        return 6 * (2 * B.radius)^2
    end
end
