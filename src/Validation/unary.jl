function validate_area(X::LazySet)
    return validate_dim2(X, 2; fun=area)
end
