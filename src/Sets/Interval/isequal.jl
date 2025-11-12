function ==(X::Interval, Y::Interval)
    return IA.isequal_interval(X.dat, Y.dat)
end
