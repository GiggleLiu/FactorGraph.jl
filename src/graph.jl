using TupleTools
using OMEinsum

export FGraph, legmap, nv, ne, vertices, edges, count_vertices, dangling_legs, occupied_legs, neighbors, add_vertex
export isvoid, isloop, isbond, isstar, is_directed, is_simple
export contract
# if allowing star like bond, the graph is not simple
# also, a tensor network can be a multi-graph, it can not be avoided in the contraction process
# note: not all legs form bonds
# note: single site bond means trace
# note: with star contraction, the relation between elimination ordering and contraction ordering still holds
# note: general graph means the equivalence between edges and vertices

struct NotSimpleGraphError <: Exception
    msg::String
end

struct FGraph{T, TT<:AbstractArray{T}, ET, VTT<:Union{NTuple{<:Any, TT}, Vector{TT}}}
    tensors::VTT
    legmap::Matrix{ET}
end

legmap(gtn::FGraph) = gtn.legmap

deleteat(a::AbstractVector, indices::Int) = deleteat!(copy(a), indices)
deleteat(a::AbstractVector, indices::Vector) = deleteat!(copy(a), indices |> sort)
deleteat(a::AbstractVector, indices::Tuple) = deleteat!(copy(a), indices |> TupleTools.sort)
function deleteat(a::AbstractMatrix, axis::Int, indices)
    Ni = length(indices)
    sizes = [size(a)...]
    sizes[axis] -= Ni
    b = similar(a, sizes...)
    inds = setdiff(1:size(a, axis), indices)
    if axis == 1
        for (i, k) in enumerate(inds)
            @inbounds b[i,:] = view(a,k,:)
        end
    else
        for (i, k) in enumerate(inds)
            @inbounds b[:,i] = view(a,:,k)
        end
    end
    b
end

insert(a::AbstractVector, i, elem) = insert!(copy(a), i, elem)

nv(gtn::FGraph) = size(gtn.legmap, 1)
ne(gtn::FGraph) = size(gtn.legmap, 2)
vertices(gtn::FGraph) = 1:nv(gtn)
edges(gtn::FGraph) = 1:ne(gtn)

# performance is not good
vertices(gtn::FGraph, ie::Int) = findall(!iszero, @inbounds view(gtn |> legmap, :, ie))
edges(gtn::FGraph, it::Int) = findall(!iszero, @inbounds view(gtn |> legmap, it, :))

isvoid(gtn::FGraph, ie::Int) = count_vertices(gtn, ie) == 0
isloop(gtn::FGraph, ie::Int) = count_vertices(gtn, ie) == 1
isbond(gtn::FGraph, ie::Int) = count_vertices(gtn, ie) == 2
isstar(gtn::FGraph, ie::Int) = count_vertices(gtn, ie) > 2

count_vertices(gtn::FGraph, ie::Int) = count(x->x!=0, view(gtn |> legmap, :,ie))

is_directed(gtn::FGraph) = false

function neighbors(gtn::FGraph, it)
    egs = vcat(edges.(Ref(gtn), it)...)
    nbs = Int[]
    for iedge in egs
        for jt in vertices(gtn, iedge)
            if !(jt in it) && !(jt in nbs)
                push!(nbs, jt)
            end
        end
    end
    nbs
end

# remove edge i
function rem_edge(gtn::FGraph, ie)
    FGraph(gtn.tensors |> copy, deleteat(gtn.legmap, 2, ie))
end

function rem_vertex(gtn::FGraph, it)
    FGraph(deleteat(gtn.tensors, it), deleteat(gtn.legmap, 1, it))
end

function rem_leg!(gtn::FGraph, iv::Int, idim::Int)
    @inbounds for ie in 1:ne(gtn)
        if gtn.legmap[iv, ie] == idim
            gtn.legmap[iv, ie] = 0
        elseif gtn.legmap[iv,ie] > idim
            gtn.legmap[iv, ie] -= 1
        end
    end
    gtn
end

function add_vertex(gtn::FGraph{T, TT}, ts::TT, edges) where {T, TT<:AbstractArray{T}}
    lm = similar(gtn |> legmap, nv(gtn)+1, ne(gtn))
    @inbounds copy!(view(lm, 1:nv(gtn), :), gtn |> legmap)
    lm[end, :] .= 0
    for (i, ie) in enumerate(edges)
        @inbounds ie isa Nothing || (lm[end, ie] = i)
    end
    FGraph(push!(copy(gtn.tensors), ts), lm)
end

for FUNC in [:rem_edge, :rem_vertex, :add_vertex, :contract]
    @eval $FUNC(args...) = gtn::FGraph -> $FUNC(gtn, args...)
    #@eval $FUNC(it::Union{Int, Vector, }) = gtn -> $FUNC(gtn, it)
end

function Base.show(io::IO, g::FGraph{T}) where T
    dir = is_directed(g) ? "directed" : "undirected"
    print(io, "{$(nv(g)), $(ne(g))} $dir $T tensor network")
    for ie in edges(g)
        print(io, "\nE($ie) = $(join(vertices(g, ie), "-"))")
    end
end
Base.show(io::IO, ::MIME"text/plain", g::FGraph) = show(io, g)

### INTERFACE
Base.eltype(x::FGraph{T}) where T = T
is_simple(gtn::FGraph) = all(ie->isbond(gtn, ie), gtn |> edges)

function occupied_legs(gtn::FGraph, it::Int)
    legmap(gtn)[it, edges(gtn, it)]
end

function dangling_legs(gtn::FGraph, it::Int)
    setdiff(1:ndims(gtn.tensors[it]), occupied_legs(gtn, it))
end

for FUNC in [:dangling_legs, :occupied_legs]
    @eval $FUNC(gtn::FGraph) = $FUNC.(Ref(gtn), gtn |> vertices)
end

for FUNC in [:count_vertices]
    @eval $FUNC(gtn::FGraph) = $FUNC.(Ref(gtn), gtn |> edges)
end

"""
contract a bond.
"""
function contract(gtn::FGraph, ie::Int)
    Ne = ne(gtn)
    vs = ie isa Tuple ? union(vertices.(Ref(gtn), ie)...) : vertices(gtn, ie)
    Nvs = length(vs)
    Nvs == 0 && return gtn |> rem_edge(ie)
    Nvs == 1 && return _contract_loop(gtn, ie, vs[])
    #Nvs != 2 && throw(NotSimpleGraphError("Star contraction is not implemented yet: $gtn"))
    IVS = edges.(Ref(gtn), vs)
    #ICON = intersect(IV...)
    IALL = union(IVS...)
    IC = setdiff(IALL, ie)
    ID = []
    #IDS = dangling_legs.(Ref(gtn), vs)
    k = Ne
    for i = 1:Nvs
        IV, v = IVS[i], vs[i]
        Nd = ndims(gtn.tensors[v])
        _code = zeros(Int, Nd)
        for ie in IV
            @inbounds _code[legmap(gtn)[v,ie]] = ie
        end
         @inbounds for idim = 1:Nd
            if _code[idim] == 0
                k += 1
               _code[idim] = k
               push!(ID, k)
            end
        end
        IVS[i] = _code
    end

    #new_tensor = eincontract((i%2==1 ? gtn.tensors[vs[(i+1)÷2]] : Tuple(IVS[(i+1)÷2]) for i=1:2*length(IVS))..., union(IC, ID)|>Tuple)
    xs = Tuple(gtn.tensors[vs])
    ixs = Tuple(Tuple.(IVS))
    iy = Tuple(union(IC, ID))
    new_tensor = einsum(ixs, xs, iy)
    NIC = Tuple(i > ie ? i-1 : i for i in IC)
    gtn |> rem_edge(ie) |> rem_vertex(vs) |> add_vertex(new_tensor, NIC)
end

function _contract_loop(gtn::FGraph, ie::Int, iv::Int)
    il = legmap(gtn)[iv, ie]
    gtn = gtn |> rem_edge(ie)
    gtn.tensors[iv] = dropdims(sum(gtn.tensors[iv], dims=il), dims=il)
    rem_leg!(gtn, iv, il)
end


# check tensor shapes, leg indices
function check_tensors(gtn::FGraph)
end

# TODO line graph
