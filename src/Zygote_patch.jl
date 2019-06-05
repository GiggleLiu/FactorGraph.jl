using Zygote: @adjoint, gradient, @adjoint!
using OMEinsum: outputtensor

@adjoint setdiff(args...) = setdiff(args...), _ -> nothing
@adjoint intersect(args...) = intersect(args...), _ -> nothing

nograd(f, args...) = f(args...)
@adjoint nograd(f, args...) = f(args...), _ -> nothing
@adjoint optimaltree(args...) = optimaltree(args...), _ -> nothing

@adjoint get_size_dict(args...) = get_size_dict(args...), _ -> nothing
@adjoint leg_analysis(args...) = leg_analysis(args...), _ -> nothing

@adjoint! function copyto!(xs::AbstractArray, ys::AbstractArray)
  xs_ = copy(xs)
  copyto!(xs, ys), function (dxs)
    copyto!(xs_, xs)
    (nothing, dxs)
  end
end
