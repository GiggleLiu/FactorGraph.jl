export asarray
asarray(x::Number) = fill(x, ())
asarray(x::AbstractArray) = x

deleteat(a::AbstractVector, indices::Int) = deleteat!(copy(a), indices)
deleteat(a::AbstractVector, indices::Vector) = deleteat!(copy(a), indices |> sort)
deleteat(a::AbstractVector, indices::Tuple) = deleteat!(copy(a), indices |> TupleTools.sort)
function deleteat(a::AbstractMatrix, axis::Int, indices)
    Ni = length(indices)
    sizes = [size(a)...]
    sizes[axis] -= Ni
    b = similar(a, sizes...)
    inds = setdiff(1:size(a, axis), indices)
    if axis == 1
        for (i, k) in enumerate(inds)
            @inbounds b[i,:] = view(a,k,:)
        end
    else
        for (i, k) in enumerate(inds)
            @inbounds b[:,i] = view(a,:,k)
        end
    end
    b
end

insert(a::AbstractVector, i, elem) = insert!(copy(a), i, elem)
