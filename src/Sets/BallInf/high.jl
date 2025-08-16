function high(B::BallInf)
    return _high_AbstractBallp(B)
end

@validate function high(B::BallInf, i::Int)
    return _high_AbstractBallp(B, i)
end
