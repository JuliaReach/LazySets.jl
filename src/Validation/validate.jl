function _validate(ex::Expr)
    # get function expression
    def = ExprTools.splitdef(ex)
    function_name = def[:name]
    function_args = def[:args]
    body_args = def[:body].args

    # get validation function and arguments to pass
    validation_fun, arg_indices = VALIDATE_DICT[function_name]

    # create assertion call of validation function
    arguments = [_unpack_arg(function_args, arg) for arg in arg_indices]
    assertion = :(@assert $validation_fun($(arguments...)))

    # set source-code annotation to first line (only useful for debugging)
    assertion.args[2] = body_args[1]

    # insert assertion in source code
    insert!(body_args, 2, assertion)

    # create new function expression
    return def
end

macro validate(ex)
    def = _validate(ex)
    f2 = ExprTools.combinedef(def)
    return quote
        $(esc(f2))
    end
end

macro validate_commutative(ex)
    # apply @validate effect
    def = _validate(ex)

    # apply @commutative effect
    f2 = ReachabilityBase.Commutative.commutative(def)

    return quote
        $(f2)
    end
end

function _unpack_arg(args::Vector, idx::Int)
    return _unpack_arg(args[idx])
end

function _unpack_arg(args::Vector, param::Symbol)
    return param
end

function _unpack_arg(arg::Symbol)
    return arg
end

function _unpack_arg(arg::Expr)
    if arg.head == :(::)
        # var::Type
        return arg.args[1]
    elseif arg.head == :kw
        # var::Type=value
        return _unpack_arg(arg.args[1])
    end
    throw(ArgumentError("unsupported argument $arg in validation"))
end
