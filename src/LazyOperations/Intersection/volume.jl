function volume(cap::Intersection)
    return volume(intersection(cap.X, cap.Y))
end
