function area(B::BallInf)
    @assert dim(B) == 2 "this function only applies to two-dimensional sets, " *
                        "but the given set is $(dim(B))-dimensional"
    return (2 * B.radius)^2
end
