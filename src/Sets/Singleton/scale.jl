function scale!(α::Real, S::Singleton)
    S.element .*= α
    return S
end
