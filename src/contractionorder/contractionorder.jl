module ContractionOrder

export ContractionTree

struct ContractionTree
    left
    right
end

include("incidencelist.jl")
include("greedy.jl")

end
