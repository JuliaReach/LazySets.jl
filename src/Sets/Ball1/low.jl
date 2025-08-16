function low(B::Ball1)
    return _low_AbstractBallp(B)
end

@validate function low(B::Ball1, i::Int)
    return _low_AbstractBallp(B, i)
end
