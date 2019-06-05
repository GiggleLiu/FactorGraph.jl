module FactorGraph

using TupleTools
using OMEinsum
using TensorOperations: optimaltree, Power
using TensorOperations

include("utils.jl")
include("contract.jl")
include("incidence_matrix.jl")
include("fgraph.jl")
include("random_graphs.jl")
include("Zygote_patch.jl")

end # module
