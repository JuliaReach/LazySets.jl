using LazySets, Test
import Aqua

@testset "Aqua tests" begin
    Aqua.test_all(LazySets; ambiguities=false,
                  # known piracies
                  piracies=(treat_as_own=[rand],))

    # do not warn about ambiguities in dependencies
    Aqua.test_ambiguities(LazySets)
end
