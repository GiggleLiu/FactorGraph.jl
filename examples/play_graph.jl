using FactorGraph, Random

# generate a random with 5 nodes (factors) and 10 edges (variables)
Random.seed!(2); rg = random_fg(5, 10)

# remove the first two edges from this incidence list, notice these labels (numbers) of edges are dynamicaly assigned.
rg = rg |> rem_edge(1) |> rem_edge(1)
# remove the first vertex from this incidence list
rg = rg |> rem_vertex(2)

# if we want to eliminate a variable (sum over an edge)
rg = rg |> contract(3) |> contract(1)
