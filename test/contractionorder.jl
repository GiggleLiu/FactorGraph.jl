import FactorGraph
using FactorGraph.ContractionOrder
using FactorGraph.ContractionOrder: analyze_contraction, contract_pair!, evaluate_costs
using OMEinsum

using Test
@testset "analyze contraction" begin
    incidence_list = IncidenceList(Dict('A' => ['a', 'b', 'k', 'o', 'f'], 'B'=>['a', 'c', 'd', 'm', 'f'], 'C'=>['b', 'c', 'e', 'f'], 'D'=>['e'], 'E'=>['d', 'f']), openedges=['c', 'f', 'o'])
    info = analyze_contraction(incidence_list, 'A', 'B')
    @test Set(info.l1) == Set(['k'])
    @test Set(info.l2) == Set(['m'])
    @test Set(info.l12) == Set(['a'])
    @test Set(info.l01) == Set(['b','o'])
    @test Set(info.l02) == Set(['c', 'd'])
    @test Set(info.l012) == Set(['f'])
end

@testset "tree greedy" begin
    incidence_list = IncidenceList(Dict('A' => ['a', 'b'], 'B'=>['a', 'c', 'd'], 'C'=>['b', 'c', 'e', 'f'], 'D'=>['e'], 'E'=>['d', 'f']))
    log2_edge_sizes = Dict([c=>i for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)
    il = copy(incidence_list)
    contract_pair!(il, 'A', 'B', log2_edge_sizes)
    target = IncidenceList(Dict('A' => ['b', 'c', 'd'], 'C'=>['b', 'c', 'e', 'f'], 'D'=>['e'], 'E'=>['d', 'f']))
    @test il.v2e == target.v2e
    @test length(target.e2v) == length(il.e2v)
    for (k,v) in il.e2v
        @test sort(target.e2v[k]) == sort(v)
    end
    costs = evaluate_costs(MinSpaceOut(), incidence_list, log2_edge_sizes)
    @test costs == Dict(('A', 'B')=>9-0.01, ('A', 'C')=>15-0.02, ('B','C')=>18-0.03, ('B','E')=>10-0.04, ('C','D')=>11-0.05, ('C', 'E')=>14-0.06)
    tree, log2_tcs, log2_scs = tree_greedy(incidence_list, log2_edge_sizes)
    @test log2_tcs == [10.0, 16, 15, 9]
    @test tree == ContractionTree(ContractionTree('A', 'B'), ContractionTree(ContractionTree('C', 'D'), 'E'))
    @test timespace_complexity(incidence_list, tree, log2_edge_sizes) == (log2(exp2(10)+exp2(16)+exp2(15)+exp2(9)), 11)
    vertices = ['A', 'B', 'C', 'D', 'E']
    optcode1 = FactorGraph.parse_eincode(incidence_list, tree, vertices=vertices)
    @test optcode1 isa OMEinsum.NestedEinsum
    tree2 = FactorGraph.parse_tree(optcode1, vertices)
    @test tree2 == tree

    eincode = ein"ab,acd,bcef,e,df->"
    size_dict = Dict([c=>(1<<i) for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)
    optcode2 = FactorGraph.optimizecode_greedy(eincode, size_dict) 
    tc, sc = timespace_complexity(optcode2, log2_edge_sizes)
    @test tc ≈ log2(exp2(10)+exp2(16)+exp2(15)+exp2(9))
    @test sc == 11
    @test optcode1 == optcode2
end