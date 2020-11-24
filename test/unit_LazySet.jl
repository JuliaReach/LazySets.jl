for ST in LazySets.subtypes(LazySet, true)
    # check that certain methods are implemented
    isconvextype(ST)
    isoperationtype(ST)
end
