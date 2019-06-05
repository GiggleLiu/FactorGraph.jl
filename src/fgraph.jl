export FGraph, legmap
export dangling_legs, occupied_legs, rem_leg!
export eliminate
export line_graph

mutable struct FGraph{T, TT<:AbstractArray{T}, ET, VTT<:Union{NTuple{<:Any, TT}, Vector{TT}}} <: IncidenceMatrix
    tensors::VTT
    legmap::Matrix{ET}
end

legmap(fg::FGraph) = fg.legmap
Base.eltype(x::FGraph{T}) where T = T
Base.copy(fg::FGraph) = FGraph(fg.tensors[:], fg.legmap)

################# IncidenceMatrix interfaces #################
imap(fg::FGraph) = legmap(fg)

# remove edge i
function rem_edge!(fg::FGraph, ie)
    fg.legmap = deleteat(fg.legmap, 2, ie)
    fg
end

function rem_vertex!(fg::FGraph, it)
    deleteat!(fg.tensors, it)
    fg.legmap = deleteat(fg.legmap, 1, it)
    fg
end

function add_vertex!(fg::FGraph{T, TT}, ts::TT, edges) where {T, TT<:AbstractArray{T}}
    lm = similar(fg |> legmap, nv(fg)+1, ne(fg))
    @inbounds copy!(view(lm, 1:nv(fg), :), fg |> legmap)
    lm[end, :] .= 0
    for (i, ie) in enumerate(edges)
        @inbounds ie isa Nothing || (lm[end, ie] = i)
    end
    fg.legmap = lm
    push!(fg.tensors, ts)
    fg
end

### NEW INTERFACE
"""remove a leg"""
function rem_leg!(fg::FGraph, iv::Int, idim::Int)
    @inbounds for ie in 1:ne(fg)
        if fg.legmap[iv, ie] == idim
            fg.legmap[iv, ie] = 0
        elseif fg.legmap[iv,ie] > idim
            fg.legmap[iv, ie] -= 1
        end
    end
    fg
end

function occupied_legs(fg::FGraph, it::Int)
    legmap(fg)[it, edges(fg, it)]
end

function dangling_legs(fg::FGraph, it::Int)
    setdiff(1:ndims(fg.tensors[it]), occupied_legs(fg, it))
end

for FUNC in [:dangling_legs, :occupied_legs]
    @eval $FUNC(fg::FGraph) = $FUNC.(Ref(fg), fg |> vertices)
end

"""
    eliminate([fg::FGraph], ie::Int) -> FGraph

eliminate a variable.
"""
function eliminate(fg::FGraph, ie::Int)
    Ne = ne(fg)
    vs = ie isa Tuple ? union(vertices.(Ref(fg), ie)...) : vertices(fg, ie)
    Nvs = length(vs)
    Nvs == 0 && return fg |> rem_edge(ie)
    Nvs == 1 && return _eliminate_loop(fg, ie, vs[])
    #Nvs != 2 && throw(NotSimpleGraphError("Star contraction is not implemented yet: $fg"))
    IVS = edges.(Ref(fg), vs)
    #ICON = intersect(IV...)
    IALL = union(IVS...)
    IC = setdiff(IALL, ie)
    ID = []
    #IDS = dangling_legs.(Ref(fg), vs)
    k = Ne
    for i = 1:Nvs
        IV, v = IVS[i], vs[i]
        Nd = ndims(fg.tensors[v])
        _code = zeros(Int, Nd)
        for ie in IV
            @inbounds _code[legmap(fg)[v,ie]] = ie
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

    #new_tensor = eincontract((i%2==1 ? fg.tensors[vs[(i+1)÷2]] : Tuple(IVS[(i+1)÷2]) for i=1:2*length(IVS))..., union(IC, ID)|>Tuple)
    xs = Tuple(fg.tensors[vs])
    ixs = Tuple(Tuple.(IVS))
    iy = Tuple(union(IC, ID))
    new_tensor = einsum(ixs, xs, iy)
    NIC = Tuple(i > ie ? i-1 : i for i in IC)
    fg |> rem_edge(ie) |> rem_vertex(vs) |> add_vertex(new_tensor, NIC)
end
eliminate(args...) = fg::FGraph -> eliminate(fg, args...)

function _eliminate_loop(fg::FGraph, ie::Int, iv::Int)
    il = legmap(fg)[iv, ie]
    fg = fg |> rem_edge(ie)
    fg.tensors[iv] = dropdims(sum(fg.tensors[iv], dims=il), dims=il)
    rem_leg!(fg, iv, il)
end

# check tensor shapes, leg indices
function check_tensors(fg::FGraph)
end

# TODO line graph
function line_graph(fg::FGraph, edge_tensors)
    FGraph(edge_tensors, fg.legmap')
end
