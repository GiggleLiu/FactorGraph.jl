export optcontract, treecontract

"""
find the optimal contraction tree.
"""
function TensorOperations.optimaltree(network, size_dict::Dict{Int, Int}=Dict{Int, Int}())
    unique_tokens = union(network...)
    optimaltree(network, Dict(token=>Power{:χ}(get(size_dict, token, 1),1) for token in unique_tokens))
end

function leg_analysis(IVS...)
    IALL = union(IVS...)
    II = intersect(IVS...)
    IC = setdiff(IALL, II)
    IALL, II, IC
end

_treecontract(tree::Int, ixs, xs, iy::Nothing) = xs[tree], ixs[tree]
function _treecontract(tree::Int, ixs, xs, iy)
    iy0, y = ixs[tree], xs[tree]
    einsum((iy0,), (C,), iy), iy
end

function _treecontract(tree, ixs, xs, iy)
    i, j = tree
    A, IA = _treecontract(i, ixs, xs, nothing)
    B, IB = _treecontract(j, ixs, xs, nothing)
    _iy = iy == nothing ? Tuple(leg_analysis(IA, IB)[3]) : iy
    einsum((IA, IB), (A, B), _iy), _iy
end

function treecontract(tree, ixs, xs, iy)
    _treecontract(tree, ixs, xs, iy) |> first
end

function get_size_dict(ixs, xs)
    nt = length(ixs)
    size_dict = Dict{Int, Int}()
    @inbounds for i = 1:nt
        for (N, leg) in zip(size(xs[i]), ixs[i])
            if haskey(size_dict, leg)
                size_dict[leg] == N || throw(DimensionMismatch("size of contraction leg $leg not match."))
            else
                size_dict[leg] = N
            end
        end
    end
    return size_dict
end

function optcontract(ixs, xs, iy)
    size_dict = get_size_dict(ixs, xs)

    print("finding optimal tree ...")
    tree, cost = optimaltree(ixs, size_dict)
    @show tree, cost
    treecontract(tree, ixs, xs, iy)
end
