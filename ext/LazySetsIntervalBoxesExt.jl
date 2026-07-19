module LazySetsIntervalBoxesExt

using IntervalArithmetic: interval
using IntervalBoxes: IntervalBox
using LazySets: AbstractHyperrectangle, UnionSetArray, dim, high, low
using LazySets.EmptySetModule: EmptySet
using LazySets.HyperrectangleModule: Hyperrectangle
import Base: convert
import LazySets: _difference

include("IntervalBoxes/HyperrectangleExt.jl")

function convert(::Type{IntervalBox}, H::AbstractHyperrectangle)
    return IntervalBox(interval.(low(H), high(H))...)
end

function _difference(X::AbstractHyperrectangle, Y::AbstractHyperrectangle)
    Xib = convert(IntervalBox, X)
    Yib = convert(IntervalBox, Y)
    Zibs = setdiff(Xib, Yib)
    if isempty(Zibs)
        return EmptySet(dim(X))
    end
    return UnionSetArray(convert.(Hyperrectangle, Zibs))
end

end  # module
