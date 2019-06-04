export sequential_contract

function _contraction_code(fg::FGraph, i::Int, idn::Int)
    DIM = ndims(fg.tensors[i])
    ids = zeros(Int, DIM)
    for ie=1:ne(fg)
        leg = fg.legmap[i,ie]
        if leg != 0
            @inbounds ids[leg] = ie
        end
    end
    for i=1:DIM
        if ids[i] == 0
            idn += 1
            ids[i] = idn
        end
    end
    ids, idn
end

"""get the contraction code from a factor graph"""
function contraction_code(fg::FGraph)
    codes = Any[]
    idn = ne(fg)
    for i in 1:nv(fg)
        code, idn = _contraction_code(fg, i, idn)
        push!(codes, code)
    end
    IC = FactorGraph.leg_analysis(codes)[3]
    codes, IC
end

using TensorOperations: BinaryTreeNode
order2tree(order) = length(order) == 1 ? order[1] : BinaryTreeNode(order2tree(order[1:end-1]), order[end])

@adjoint contraction_code(fg) = contraction_code(fg), adjy->nothing
@adjoint order2tree(order) = order2tree(order), adjy->nothing

"""
    sequential_contract(fg::FGraph, order)

Contract tensors in a factor graph by `order` of tensors.
"""
function sequential_contract(fg::FGraph, order)
    ixs, iy = contraction_code(fg)
    tss = Tuple(fg.tensors)
    tree = order2tree(order)
    treecontract(tree, Tuple.(ixs), tss, Tuple(iy))
end
