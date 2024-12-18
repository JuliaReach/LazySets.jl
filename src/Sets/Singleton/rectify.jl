function rectify(S::Singleton)
    return Singleton(rectify(element(S)))
end
