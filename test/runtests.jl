using FactorGraph
using Test

@testset "graph" begin
    include("graph.jl")
end

@testset "contract" begin
    include("contract.jl")
end
