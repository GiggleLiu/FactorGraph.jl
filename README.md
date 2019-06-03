# FactorGraph

Incidence list base factor graph modeling.

WIP...

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GiggleLiu.github.io/FactorGraph.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GiggleLiu.github.io/FactorGraph.jl/dev)
[![Build Status](https://travis-ci.com/GiggleLiu/FactorGraph.jl.svg?branch=master)](https://travis-ci.com/GiggleLiu/FactorGraph.jl)
[![Codecov](https://codecov.io/gh/GiggleLiu/FactorGraph.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/GiggleLiu/FactorGraph.jl)

```julia console
julia> using FactorGraph, Random

julia> # generate a random with 5 nodes (factors) and 10 edges (variables)
       Random.seed!(2); rg = random_fg(5, 10)
FGraph{Float64}(V5,E10) f1(1,2,3,4,5) * f2(3) * f3(1,2,5) * f4(3,4) * f5(3) * f6(2,4,5) * f7(3,4) * f8(3) * f9() * f10(2,3,4)

julia> # remove the first two edges from this incidence list, notice these labels (numbers) of edges are dynamicaly assigned.
       rg = rg |> rem_edge(1) |> rem_edge(1)
FGraph{Float64}(V5,E8) f1(1,2,5) * f2(3,4) * f3(3) * f4(2,4,5) * f5(3,4) * f6(3) * f7() * f8(2,3,4)

julia> # remove the first vertex from this incidence list
       rg = rg |> rem_vertex(2)
FGraph{Float64}(V4,E8) f1(1,4) * f2(2,3) * f3(2) * f4(3,4) * f5(2,3) * f6(2) * f7() * f8(2,3)

julia> # eliminate a variable (sum over an edge)
       rg = rg |> contract(3) |> contract(1)
FGraph{Float64}(V3,E6) f1(1,2) * f2(2,3) * f3(1,2) * f4(1) * f5() * f6(1,2)

julia> # the resulting tensors
       rg.tensors .|> size
3-element Array{Tuple{Int64,Int64,Int64,Vararg{Int64,N} where N},1}:
 (2, 2, 2, 2, 2, 2)
 (2, 2, 2, 2, 2)   
 (2, 2, 2)         

julia> # the topology of the incidence list of the resulting V-3, E-4 graph.
       rg.legmap
3×4 Array{Int64,2}:
 4  5  0  6
 4  0  0  5
 0  0  0  0
```
