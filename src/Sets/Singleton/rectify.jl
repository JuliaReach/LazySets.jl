function rectify(S::Singleton)
    return Singleton(rectify(center(S)))
end
