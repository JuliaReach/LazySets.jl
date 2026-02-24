@validate function area(B::Ball2)
    if dim(B) == 2
        return π * B.radius^2
    else
        return 4 * π * B.radius^2
    end
end
