function scale(α::Real, S::Singleton)
    return _scale_copy_inplace(α, S)
end

function scale!(α::Real, S::Singleton)
    S.element .*= α
    return S
end
