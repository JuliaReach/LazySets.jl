# since there is only a single instance (per dimensionality), we do not copy
function copy(Z::ZeroSet)
    return Z
end
