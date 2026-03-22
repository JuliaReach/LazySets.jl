function copy(S::Singleton)
    return Singleton(copy(S.element))
end
