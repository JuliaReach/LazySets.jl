for ST in LazySets.subtypes(ConvexSet, true)
    # check that certain methods are implemented
    isconvextype(ST)
    isoperationtype(ST)
    isboundedtype(ST)
end
