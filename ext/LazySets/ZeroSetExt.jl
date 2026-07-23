using LazySets.SingletonModule: Singleton
using LazySets.ZeroSetModule: ZeroSet
import LazySets.API: translate

@validate function translate(Z::ZeroSet, v::AbstractVector)
    return Singleton(v)
end
