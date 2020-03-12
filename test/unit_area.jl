# not implemented for dimension != 2
p = BallInf(zeros(3), 1.0)
@test_throws ArgumentError surface(p)
@test_throws AssertionError area(p)

# sets with zero area
p = Singleton([0.0, 1.0])
@test area(p) ≈ 0.0
p = Interval(0.0, 1.0) × Interval(0.0, 0.0)
@test area(p) ≈ 0.0

# triangle
p = VPolygon([[13, 14.], [16, 30.], [50, 10.]])
@test area(p) ≈ 302.0

# quadrilateral
vlist = [[3, 4.], [5, 11.], [12, 8.], [9, 5.], [5, 6.]]
p = VPolygon(vlist) # convex hull => 4 points
@test area(p) ≈ 35.0

# test surface function, which falls back to the area function for 2D sets
@test surface(p) ≈ 35.0

# general case (non-convex)
@test LazySets._area_polygon(vlist) ≈ 30.0
