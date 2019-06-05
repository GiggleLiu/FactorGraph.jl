#=
NOTE:
if allowing star like bond, the graph is not simple
also, a tensor network can be a multi-graph, it can not be avoided in the contraction process
note: not all legs form bonds
note: single site bond means trace
note: with star contraction, the relation between elimination ordering and contraction ordering still holds
note: general graph means the equivalence between edges and vertices
=#

export IncidenceMatrix, imap
export linegraph_imap
export nv, ne, vertices, edges, count_vertices, neighbors, findedges
export add_vertex, rem_edge, rem_vertex, add_vertex!, rem_edge!, rem_vertex!
export isvoid, isloop, isbond, isstar, isdirected_graph, issimple_graph, isconnected
export NotSimpleGraphError

abstract type IncidenceMatrix end

########### interfaces ##############
"""
    imap(g::IncidenceMatrix) -> AbstractArray

get the incidence matrix"""
function imap end

"""remove an edge"""
function rem_edge! end

"""remove a vertex"""
function rem_vertex! end

"""add a new vertex"""
function add_vertex! end

Base.show(io::IO, g::IncidenceMatrix) = show(io, "plain/text", g)
function Base.show(io::IO, ::MIME"plain/text", g::IncidenceMatrix)
    print(io, "$(typeof(g).name)")
    print(io, "(V$(nv(g)),E$(ne(g))) ")
    print(io, join(["f$ie($(join(vertices(g, ie), ",")))" for ie in edges(g)], " * "))
end

# NOTE: `copy` in base is also required.

############ derived fallback functions with above interfaces ##############
for FUNC in [:rem_edge, :rem_vertex, :add_vertex]
    IFUNC = Symbol(FUNC, :!)
    @eval $FUNC(fg::IncidenceMatrix, args...) = $IFUNC(copy(fg), args...)
    @eval $FUNC(args...) = fg::IncidenceMatrix -> $FUNC(fg, args...)
    @eval $IFUNC(args...) = fg::IncidenceMatrix -> $IFUNC(fg, args...)
end

nv(fg::IncidenceMatrix) = size(fg |> imap, 1)
ne(fg::IncidenceMatrix) = size(fg |> imap, 2)

vertices(fg::IncidenceMatrix) = 1:nv(fg)
vertices(fg::IncidenceMatrix, ie::Int) = findall(!iszero, @inbounds view(fg |> imap, :, ie))
edges(fg::IncidenceMatrix) = 1:ne(fg)
edges(fg::IncidenceMatrix, it::Int) = findall(!iszero, @inbounds view(fg |> imap, it, :))

count_vertices(fg::IncidenceMatrix, ie::Int) = count(x->x!=0, view(fg |> imap, :,ie))
for FUNC in [:count_vertices]
    @eval $FUNC(fg::IncidenceMatrix) = $FUNC.(Ref(fg), fg |> edges)
end

# boolean functions about local features
isvoid(fg::IncidenceMatrix, ie::Int) = count_vertices(fg, ie) == 0
isloop(fg::IncidenceMatrix, ie::Int) = count_vertices(fg, ie) == 1
isbond(fg::IncidenceMatrix, ie::Int) = count_vertices(fg, ie) == 2
isstar(fg::IncidenceMatrix, ie::Int) = count_vertices(fg, ie) > 2
"""
    isconnected(fg::IncidenceMatrix, vertices, [edge]) -> Bool

Return true if vertices are connected (by edge).
"""
isconnected(fg::IncidenceMatrix, vertices, edge) = all(vertex->imap(fg)[vertex, edge] .!= 0, vertices)
isconnected(fg::IncidenceMatrix, vertices) = any(edge->isconnected(fg, vertices, edge), edges(fg))

# boolean function about global features
isdirected_graph(fg::IncidenceMatrix) = false
issimple_graph(fg::IncidenceMatrix) = all(ie->isbond(fg, ie), fg |> edges)

"""
    neighbors(fg::IncidenceMatrix, it) -> Vector

the neighbors of a vertex.
"""
function neighbors(fg::IncidenceMatrix, it)
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

linegraph_imap(fg::IncidenceMatrix) = imap(fg)'

"""
    findedges([func], fg::IncidenceMatrix, vertices) -> edge(s)

find edges relating vertices, alternative argument `func` is the finder, it can be `findfirst`, `findall`...
"""
findedges(func, fg::IncidenceMatrix, vertices) = func(edge->isconnected(fg, vertices, edge), edges(fg))
findedges(fg::IncidenceMatrix, vertices) = findedges(findall, fg, vertices)

######### error handling and warning  ############
struct NotSimpleGraphError <: Exception
    msg::String
end
