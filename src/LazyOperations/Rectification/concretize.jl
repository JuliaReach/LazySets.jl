function concretize(R::Rectification)
    return concretize(rectify(concretize(R.X)))
end
