@validate function area(B::BallInf)
    if dim(B) == 2
        return (2 * B.radius)^2
    else
        return 6 * (2 * B.radius)^2
    end
end
