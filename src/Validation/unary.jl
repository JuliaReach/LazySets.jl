function validate_area(X::LazySet)
    return validate_dims(X, (2, 3); fun=area)
end
