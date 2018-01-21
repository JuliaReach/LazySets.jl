# check for method ambiguities
@test check_method_ambiguity_binary(⊆)
@test check_method_ambiguity_binary(is_intersection_empty)
