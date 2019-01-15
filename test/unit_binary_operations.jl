# check for method ambiguities
@test check_method_ambiguity_binary(intersection)
@test check_method_ambiguity_binary(issubset)
@test check_method_ambiguity_binary(is_intersection_empty)
@test check_method_ambiguity_binary(convex_hull)
