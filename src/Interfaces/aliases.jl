export CompactSet,
       NonCompactSet

"""
    CompactSet

An alias for compact set types.

### Notes

Most lazy operations are not captured by this alias because whether their result
is compact or not depends on the argument(s).
"""
const CompactSet = Union{AbstractCentrallySymmetric,
                         AbstractPolytope,
                         AbstractPolynomialZonotope,
                         EmptySet,
                         Polygon}

"""
    NonCompactSet

An alias for non-compact set types.

### Notes

Most lazy operations are not captured by this alias because whether their result
is non-compact or not depends on the argument(s).
"""
const NonCompactSet = Union{HalfSpace,
                            Hyperplane,
                            Line,
                            Line2D,
                            Universe}
