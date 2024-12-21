function convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)
    return Hyperrectangle(center(H), radius_hyperrectangle(H))
end

function load_IntervalArithmetic_convert()
    return quote
        import IntervalArithmetic as IA

        function convert(::Type{IA.IntervalBox}, H::AbstractHyperrectangle)
            return IA.IntervalBox(IA.interval.(low(H), high(H)))
        end

        function convert(::Type{Hyperrectangle}, IB::IA.IntervalBox)
            low_IB = IA.inf.(IB)
            high_IB = IA.sup.(IB)
            return Hyperrectangle(; low=low_IB, high=high_IB)
        end

        function convert(::Type{Hyperrectangle}, I::IA.Interval)
            low_I = [IA.inf(I)]
            high_I = [IA.sup(I)]
            return Hyperrectangle(; low=low_I, high=high_I)
        end
    end
end  # quote / load_IntervalArithmetic_convert
