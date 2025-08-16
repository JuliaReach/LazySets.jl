function low(B::BallInf)
    return _low_AbstractBallp(B)
end

@validate function low(B::BallInf, i::Int)
    return _low_AbstractBallp(B, i)
end
