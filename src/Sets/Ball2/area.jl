function area(B::Ball2)
    @assert dim(B) == 2 "this function only applies to two-dimensional sets, " *
                        "but the given set is $(dim(B))-dimensional"
    return Base.pi * B.radius^2
end
