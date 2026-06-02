# faster support-function evaluation for SingleEntryVector
@validate function ρ(d::SingleEntryVector, sih::SymmetricIntervalHull)
    return abs(d.v) * get_radius!(sih, d.i)
end
