function volume(cup::UnionSet)
    return volume(cup.X) + volume(cup.Y) - volume(Intersection(cup.X, cup.Y))
end
