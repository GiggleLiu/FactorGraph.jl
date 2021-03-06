using FactorGraph
using OMEinsum
using Zygote: @adjoint
using Zygote

include("sequential.jl")

"""
convert bonds to a incidence list - legmap
"""
function legmap_frombonds(bonds)
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
@adjoint legmap_frombonds(bonds) = legmap_frombonds(bonds), adjy->nothing

"""the tensor defined on the site of a ferromagnetic ising model"""
function trg_T0(K)
    # sqrt(exp.(-K K; K -K))
    M = [sqrt(cosh(K)) sqrt(sinh(K));
         sqrt(cosh(K)) -sqrt(sinh(K))]
    T = ein"ij,ik,il->jkl"(M, M, M)
end

function fullerene(K)
    T0 = trg_T0(K)
    tensors = Tuple([T0 for i=1:60])

    # construct the factor graph
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
    lm = legmap_frombonds(bonds)
    fg = FGraph(tensors, lm)

    # get the optimal contraction tree, optmaltree does not work, unfortunately
    order = (7, 14, 13, 5, 56,
        25, 26, 58, 29, 57,
        4, 44, 36, 35, 43, 3, 21, 33, 34, 22,
        50, 38, 37, 49, 6, 31, 45, 46, 32, 8,
        47, 54, 55, 48, 12, 16, 40, 39, 15, 10,
        1, 41, 19, 17, 18, 20, 42, 2, 60, 28, 52, 30, 51, 27, 59,
        9, 23, 24, 11, 53)

    y = sequential_contract(fg, order)
    (y[] |> log)/60
end

using Test
@testset "fullerene backward" begin
    ngradient(f, x) = (f(x+1e-5) - f(x-1e-5))/2e-5
    # ferromagnetic model
    @test fullerene(1.0) ≈ 1.5147392357988083
    @test ngradient(fullerene, 0.1) ≈ fullerene'(0.1)
end
