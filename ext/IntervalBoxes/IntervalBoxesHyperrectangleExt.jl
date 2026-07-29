using IntervalArithmetic: inf, sup
using IntervalBoxes: IntervalBox
using LazySets.HyperrectangleModule: Hyperrectangle
import Base: convert

function convert(::Type{Hyperrectangle}, Ibox::IntervalBox)
    low_IB = inf.(Ibox)
    high_IB = sup.(Ibox)
    return Hyperrectangle(; low=low_IB, high=high_IB)
end
