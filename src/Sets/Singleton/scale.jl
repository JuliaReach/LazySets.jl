function scale(α::Real, S::Singleton)
    return Singleton(α * S.element)
end

function scale!(α::Real, S::Singleton)
    S.element .*= α
    return S
end
