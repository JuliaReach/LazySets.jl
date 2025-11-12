function convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)
    return Hyperrectangle(center(H), radius_hyperrectangle(H))
end

function load_IntervalArithmetic_convert()
    return quote
        import IntervalArithmetic as IA

        function convert(::Type{Hyperrectangle}, I::IA.Interval)
            low_I = [IA.inf(I)]
            high_I = [IA.sup(I)]
            return Hyperrectangle(; low=low_I, high=high_I)
        end
    end
end  # quote / load_IntervalArithmetic_convert

function load_IntervalBoxes_convert()
    return quote
        import IntervalArithmetic as IB

        function convert(::Type{IA.IntervalBox}, H::AbstractHyperrectangle)
            return IB.IntervalBox(IA.interval.(low(H), high(H)))
        end

        function convert(::Type{Hyperrectangle}, Ibox::IA.IntervalBox)
            low_IB = IB.inf.(Ibox)
            high_IB = IB.sup.(Ibox)
            return Hyperrectangle(; low=low_IB, high=high_IB)
        end
    end
end  # quote / load_IntervalBoxes_convert
