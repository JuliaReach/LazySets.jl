# if ``x = [a, b]``, then ``-x = [-b, -a]``
function reflect(X::Interval)
    return Interval(-_max(X), -_min(X))
end
