using FactorGraph.ContractionOrder
using FactorGraph.ContractionOrder: analyze_contraction, contract_pair!, evaluate_costs

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
    contract_pair!(incidence_list, 'A', 'B', log2_edge_sizes)
    target = IncidenceList(Dict('A' => ['b', 'c', 'd'], 'C'=>['b', 'c', 'e', 'f'], 'D'=>['e'], 'E'=>['d', 'f']))
    @test incidence_list.v2e == target.v2e
    @test length(target.e2v) == length(incidence_list.e2v)
    for (k,v) in incidence_list.e2v
        @test sort(target.e2v[k]) == sort(v)
    end
    incidence_list = IncidenceList(Dict('A' => ['a', 'b'], 'B'=>['a', 'c', 'd'], 'C'=>['b', 'c', 'e', 'f'], 'D'=>['e'], 'E'=>['d', 'f']))
    costs = evaluate_costs(MinSpaceOut(), incidence_list, log2_edge_sizes)
    @test costs == Dict(('A', 'B')=>9-0.01, ('A', 'C')=>15-0.02, ('B','C')=>18-0.03, ('B','E')=>10-0.04, ('C','D')=>11-0.05, ('C', 'E')=>14-0.06)
    tree, log2_tcs, log2_scs = tree_greedy(incidence_list, log2_edge_sizes)
    @test log2_tcs == [10.0, 16, 15, 9]
    @test tree == ContractionTree(ContractionTree('A', 'B'), ContractionTree(ContractionTree('C', 'D'), 'E'))
end