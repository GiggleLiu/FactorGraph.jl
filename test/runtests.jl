using FactorGraph
using Test

@testset "fgraph" begin
    include("fgraph.jl")
end

@testset "contract" begin
    include("contract.jl")
end
