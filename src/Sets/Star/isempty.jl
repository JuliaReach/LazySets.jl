function isempty(X::Star)
    return isempty(predicate(X))
end
