function convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)
    return Hyperrectangle(center(H), radius_hyperrectangle(H))
end

function load_IntervalArithmetic_convert()
    return quote
        import .IntervalArithmetic as IA

        function convert(::Type{Hyperrectangle}, I::IA.Interval)
            low_I = [IA.inf(I)]
            high_I = [IA.sup(I)]
            return Hyperrectangle(; low=low_I, high=high_I)
        end
    end
end  # quote / load_IntervalArithmetic_convert

function load_IntervalBoxes_convert_Hyperrectangle()
    return quote
        import .IntervalBoxes as IB

        function convert(::Type{Hyperrectangle}, Ibox::IB.IntervalBox)
            low_IB = IB.inf.(Ibox)  # NOTE: `inf`/`sup` are defined in IntervalArithmetic
            high_IB = IB.sup.(Ibox)
            return Hyperrectangle(; low=low_IB, high=high_IB)
        end
    end
end  # quote / load_IntervalBoxes_convert_Hyperrectangle
