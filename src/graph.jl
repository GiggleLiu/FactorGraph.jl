export FGraph, legmap, nv, ne, vertices, edges, count_vertices, dangling_legs, occupied_legs, neighbors, add_vertex, rem_edge, rem_vertex, rem_leg!
export isvoid, isloop, isbond, isstar, is_directed, is_simple
export eliminate
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

legmap(fg::FGraph) = fg.legmap

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

nv(fg::FGraph) = size(fg.legmap, 1)
ne(fg::FGraph) = size(fg.legmap, 2)
vertices(fg::FGraph) = 1:nv(fg)
edges(fg::FGraph) = 1:ne(fg)

# performance is not good
vertices(fg::FGraph, ie::Int) = findall(!iszero, @inbounds view(fg |> legmap, :, ie))
edges(fg::FGraph, it::Int) = findall(!iszero, @inbounds view(fg |> legmap, it, :))

isvoid(fg::FGraph, ie::Int) = count_vertices(fg, ie) == 0
isloop(fg::FGraph, ie::Int) = count_vertices(fg, ie) == 1
isbond(fg::FGraph, ie::Int) = count_vertices(fg, ie) == 2
isstar(fg::FGraph, ie::Int) = count_vertices(fg, ie) > 2

count_vertices(fg::FGraph, ie::Int) = count(x->x!=0, view(fg |> legmap, :,ie))

is_directed(fg::FGraph) = false

function neighbors(fg::FGraph, it)
    egs = vcat(edges.(Ref(fg), it)...)
    nbs = Int[]
    for iedge in egs
        for jt in vertices(fg, iedge)
            if !(jt in it) && !(jt in nbs)
                push!(nbs, jt)
            end
        end
    end
    nbs
end

# remove edge i
function rem_edge(fg::FGraph, ie)
    FGraph(fg.tensors |> copy, deleteat(fg.legmap, 2, ie))
end

function rem_vertex(fg::FGraph, it)
    FGraph(deleteat(fg.tensors, it), deleteat(fg.legmap, 1, it))
end

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

function add_vertex(fg::FGraph{T, TT}, ts::TT, edges) where {T, TT<:AbstractArray{T}}
    lm = similar(fg |> legmap, nv(fg)+1, ne(fg))
    @inbounds copy!(view(lm, 1:nv(fg), :), fg |> legmap)
    lm[end, :] .= 0
    for (i, ie) in enumerate(edges)
        @inbounds ie isa Nothing || (lm[end, ie] = i)
    end
    FGraph(push!(copy(fg.tensors), ts), lm)
end

for FUNC in [:rem_edge, :rem_vertex, :add_vertex, :eliminate]
    @eval $FUNC(args...) = fg::FGraph -> $FUNC(fg, args...)
    #@eval $FUNC(it::Union{Int, Vector, }) = fg -> $FUNC(fg, it)
end

### INTERFACE
Base.eltype(x::FGraph{T}) where T = T
is_simple(fg::FGraph) = all(ie->isbond(fg, ie), fg |> edges)

function occupied_legs(fg::FGraph, it::Int)
    legmap(fg)[it, edges(fg, it)]
end

function dangling_legs(fg::FGraph, it::Int)
    setdiff(1:ndims(fg.tensors[it]), occupied_legs(fg, it))
end

for FUNC in [:dangling_legs, :occupied_legs]
    @eval $FUNC(fg::FGraph) = $FUNC.(Ref(fg), fg |> vertices)
end

for FUNC in [:count_vertices]
    @eval $FUNC(fg::FGraph) = $FUNC.(Ref(fg), fg |> edges)
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

function _eliminate_loop(fg::FGraph, ie::Int, iv::Int)
    il = legmap(fg)[iv, ie]
    fg = fg |> rem_edge(ie)
    fg.tensors[iv] = dropdims(sum(fg.tensors[iv], dims=il), dims=il)
    rem_leg!(fg, iv, il)
end

Base.show(io::IO, g::FGraph) = show(io, "plain/text", g)
function Base.show(io::IO, ::MIME"plain/text", g::FGraph{T}) where T
    dir = is_directed(g) ? "directed" : "undirected"
    print(io, "FGraph{$T}")
    print(io, "(V$(nv(g)),E$(ne(g))) ")
    print(io, join(["f$ie($(join(vertices(g, ie), ",")))" for ie in edges(g)], " * "))
end

# check tensor shapes, leg indices
function check_tensors(fg::FGraph)
end

# TODO line graph
