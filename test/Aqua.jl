using LazySets, Test
import Aqua

@testset "Aqua tests" begin
    Aqua.test_all(LazySets; ambiguities=false,
                  # the unbound args should be resolved in the future
                  unbound_args=(broken=true,),
                  # known piracies
                  piracies=(treat_as_own=[rand, LazySets.activate_assertions,
                                          LazySets.deactivate_assertions],))

    # do not warn about ambiguities in dependencies
    Aqua.test_ambiguities(LazySets)
end
