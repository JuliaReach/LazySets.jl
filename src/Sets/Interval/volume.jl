function volume(X::Interval)
    return _max(X) - _min(X)
end
