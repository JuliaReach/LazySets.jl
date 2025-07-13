@validate function area(B::Ball2)
    if dim(B) == 2
        return Base.pi * B.radius^2
    else
        return 4 * Base.pi * B.radius^2
    end
end
