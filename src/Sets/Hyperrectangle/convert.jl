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
