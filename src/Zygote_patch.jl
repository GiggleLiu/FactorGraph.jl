using Zygote: @adjoint, gradient
using OMEinsum: outputtensor

@adjoint setdiff(args...) = setdiff(args...), _ -> nothing
@adjoint intersect(args...) = intersect(args...), _ -> nothing
@adjoint outputtensor(args...) = outputtensor(args...), _ -> nothing

nograd(f, args...) = f(args...)
@adjoint nograd(f, args...) = f(args...), _ -> nothing
@adjoint optimaltree(args...) = optimaltree(args...), _ -> nothing

@adjoint get_size_dict(args...) = get_size_dict(args...), _ -> nothing
@adjoint leg_analysis(args...) = leg_analysis(args...), _ -> nothing
