function high(B::Ball1)
    return _high_AbstractBallp(B)
end

@validate function high(B::Ball1, i::Int)
    return _high_AbstractBallp(B, i)
end
