function center(X::Star)
    return X.c
end

@validate function center(X::Star, i::Int)
    return X.c[i]
end
