using FactorGraph
using OMEinsum
using Zygote: @adjoint
using Zygote

bonds = ((1, 10), (1, 41), (1, 59), (2, 12), (2, 42), (2, 60), (3, 6), (3,
    43), (3, 57), (4, 8), (4, 44), (4, 58), (5, 13), (5, 56), (5,
    57), (6, 10), (6, 31), (7, 14), (7, 56), (7, 58), (8, 12), (8,
    32), (9, 23), (9, 53), (9, 59), (10, 15), (11, 24), (11, 53), (11,
    60), (12, 16), (13, 14), (13, 25), (14, 26), (15, 27), (15,
    49), (16, 28), (16, 50), (17, 18), (17, 19), (17, 54), (18,
    20), (18, 55), (19, 23), (19, 41), (20, 24), (20, 42), (21,
    31), (21, 33), (21, 57), (22, 32), (22, 34), (22, 58), (23,
    24), (25, 35), (25, 43), (26, 36), (26, 44), (27, 51), (27,
    59), (28, 52), (28, 60), (29, 33), (29, 34), (29, 56), (30,
    51), (30, 52), (30, 53), (31, 47), (32, 48), (33, 45), (34,
    46), (35, 36), (35, 37), (36, 38), (37, 39), (37, 49), (38,
    40), (38, 50), (39, 40), (39, 51), (40, 52), (41, 47), (42,
    48), (43, 49), (44, 50), (45, 46), (45, 54), (46, 55), (47,
    54), (48, 55))

function pairs2legmap(bonds)
    legmap = zeros(Int, 60, length(bonds))
    occ = zeros(Int, 60)
    for (ib,bond) in enumerate(bonds)
        A, B = bond
        occ[A] += 1
        occ[B] += 1
        legmap[A, ib] = occ[A]
        legmap[B, ib] = occ[B]
    end
    legmap
end

function trg_T0(K::Float64)
    M = [sqrt(cosh(K)) sqrt(sinh(K));
         sqrt(cosh(K)) -sqrt(sinh(K))]
    einsum(((1,2), (1,3), (1,4)), (M, M, M), (2,3,4))
end

function _code(tn::FGraph, i::Int, idn::Int)
    DIM = ndims(tn.tensors[i])
    ids = zeros(Int, DIM)
    for ie=1:ne(tn)
        leg = tn.legmap[i,ie]
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

function code(tn::FGraph)
    codes = Any[]
    idn = ne(tn)
    for i in 1:nv(tn)
        code, idn = _code(tn, i, idn)
        push!(codes, code)
    end
    IC = FactorGraph.leg_analysis(codes)[3]
    push!(codes, IC)
end

const order = (7, 14, 13, 5, 56,
    25, 26, 58, 29, 57,
    4, 44, 36, 35, 43, 3, 21, 33, 34, 22,
    50, 38, 37, 49, 6, 31, 45, 46, 32, 8,
    47, 54, 55, 48, 12, 16, 40, 39, 15, 10,
    1, 41, 19, 17, 18, 20, 42, 2, 60, 28, 52, 30, 51, 27, 59,
    9, 23, 24, 11, 53)

#optcontract(args...)
order2tree(order) = length(order) == 1 ? order[1] : [order2tree(order[1:end-1]), order[end]]
tree = order2tree(order)

using BenchmarkTools
@benchmark (treecontract(tree, args...)[] |> log)/60 seconds=1

@adjoint pairs2legmap(bonds) = pairs2legmap(bonds), adjy->nothing
@adjoint code(tn) = code(tn), adjy->(adjy,)

function fullerene(K)
    T0 = trg_T0(K)
    tensors = Tuple([T0 for i=1:60])

    lm = pairs2legmap(bonds)
    tn = FGraph(tensors, lm)

    ccode = Tuple.(code(tn))
    tss = Tuple(tn.tensors)
    Zygote.hook(adjy->(@show adjy; adjy), tss)
    y = treecontract(tree, ccode[1:end-1], tss, ccode[end])
    Zygote.hook(adjy->(@show adjy; adjy), y)
    (y[] |> log)/60
end

fullerene(1.0)
fullerene'(1.0)

using FactorGraph: random_simple_fg
using TensorOperations
using FactorGraph
g = random_simple_fg(ComplexF64, Tuple.(fill(2, 3) for i=1:30), 40, bias_factor=-5)

function network(fg::FGraph)
    [Tuple(edges(fg, ie)) for ie in vertices(fg)]
end

network(g)

net = network(g)
size_dict = FactorGraph.get_size_dict(net, g.tensors)

unique_tokens = union(net...)
c = Dict(token=>TensorOperations.Power{:χ}(get(size_dict, token, 1),1) for token in unique_tokens)
TensorOperations.optimaltree(net, c)
tree, cost = TensorOperations.optimaltree(net, size_dict)

nitem(t::TensorOperations.BinaryTreeNode) = nitem(t.left) + nitem(t.right)
nitem(x) = 1

nitem(tree)
