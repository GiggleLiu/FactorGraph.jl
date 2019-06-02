asarray(x::Number) = fill(x, ())
asarray(x::AbstractArray) = x
