module FactorGraph

using TupleTools
using OMEinsum
using TensorOperations: optimaltree, Power
using TensorOperations

include("utils.jl")
include("contract.jl")
include("graph.jl")
include("random_graphs.jl")
include("Zygote_patch.jl")

end # module
