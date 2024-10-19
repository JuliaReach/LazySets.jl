# if ``x = [a, b]``, then ``-x = [-b, -a]``
function reflect(x::Interval)
    return Interval(-max(x), -min(x))
end
