# if ``x = [a, b]``, then ``-x = [-b, -a]``
function reflect(X::Interval)
    return Interval(-max(X), -min(X))
end
