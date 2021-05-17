export tree_greedy, MinSpaceOut, MinSpaceDiff, ContractionTree

struct MinSpaceOut end
struct MinSpaceDiff end

struct LegInfo{ET}
    l1::Vector{ET}
    l2::Vector{ET}
    l12::Vector{ET}
    l01::Vector{ET}
    l02::Vector{ET}
    l012::Vector{ET}
end

struct ContractionTree
    left
    right
end

"""
    tree_greedy(incidence_list, log2_sizes)

Compute greedy order, and the time and space complexities, the rows of the `incidence_list` are vertices and columns are edges.
`log2_sizes` are defined on edges.
"""
function tree_greedy(incidence_list::IncidenceList{VT,ET}, log2_edge_sizes; method=MinSpaceOut()) where {VT,ET}
    n = nv(incidence_list)
    #tree_dict = collect(Any, 1:n)
    log2_tcs = Float64[] # time complexity
    log2_scs = Float64[]

    tree = Dict{VT,Any}([v=>v for v in vertices(incidence_list)])
    while true
        cost_values = evaluate_costs(method, incidence_list, log2_edge_sizes)
        local pair
        vmin = Inf
        for (k,v) in cost_values
            if v < vmin
                pair = k
                vmin = v
            end
        end
        log2_tc_step, code = contract_pair!(incidence_list, pair..., log2_edge_sizes)
        push!(log2_tcs, log2_tc_step)
        push!(log2_scs, space_complexity(incidence_list, log2_edge_sizes))
        if nv(incidence_list) > 1
            tree[pair[1]] = ContractionTree(tree[pair[1]], tree[pair[2]])
        else
            return ContractionTree(tree[pair[1]], tree[pair[2]]), log2_tcs, log2_scs
        end
    end
end

function contract_pair!(incidence_list, vi, vj, log2_edge_sizes)
    log2dim(legs) = sum(l->log2_edge_sizes[l], legs, init=0)
    @assert vj > vi
    # compute time complexity and output tensor
    legsets = analyze_contraction(incidence_list, vi, vj)
    D12,D01,D02,D012 = log2dim.(getfield.(Ref(legsets),3:6))
    tc = D12+D01+D02+D012  # dangling legs D1 and D2 do not contribute

    # einsum code
    eout = legsets.l01 ∪ legsets.l02 ∪ legsets.l012
    code = (edges(incidence_list, vi), edges(incidence_list, vj)) => eout

    # change incidence_list
    delete_vertex!(incidence_list, vj)
    change_edges!(incidence_list, vi, eout)
    for e in eout
        replace_vertex!(incidence_list, e, vj=>vi)
    end
    remove_edges!(incidence_list, legsets.l1 ∪ legsets.l2 ∪ legsets.l12)
    return tc, code
end

function evaluate_costs(method, incidence_list::IncidenceList{VT,ET}, log2_edge_sizes) where {VT,ET}
    # initialize cost values
    cost_values = Dict{Tuple{VT,VT},Float64}()
    for vi = vertices(incidence_list)
        for vj in neighbors(incidence_list, vi)
            if vj > vi
                cost_values[(vi,vj)] = greedy_loss(method, incidence_list, log2_edge_sizes, vi, vj)
            end
        end
    end
    return cost_values
end

function analyze_contraction(incidence_list::IncidenceList{VT,ET}, vi::VT, vj::VT) where {VT,ET}
    ei = edges(incidence_list, vi)
    ej = edges(incidence_list, vj)
    leg012,leg12,leg1,leg2,leg01,leg02 = ET[], ET[], ET[], ET[], ET[], ET[]
    # external legs
    for leg in ei ∪ ej
        isext = leg ∈ incidence_list.openedges || !all(x->x==vi || x==vj, vertices(incidence_list, leg))
        if isext
            if leg ∈ ei
                if leg ∈ ej
                    push!(leg012, leg)
                else
                    push!(leg01, leg)
                end
            else
                push!(leg02, leg)
            end
        else
            if leg ∈ ei
                if leg ∈ ej
                    push!(leg12, leg)
                else
                    push!(leg1, leg)
                end
            else
                push!(leg2, leg)
            end
        end
    end
    return LegInfo(leg1, leg2, leg12, leg01, leg02, leg012)
end

function greedy_loss(::MinSpaceOut, incidence_list, log2_edge_sizes, vi, vj)
    log2dim(legs) = sum(l->log2_edge_sizes[l], legs, init=0)
    legs = analyze_contraction(incidence_list, vi, vj)
    log2dim(legs.l01)+log2dim(legs.l02)+log2dim(legs.l012) - 0.01*(log2dim(legs.l12)+log2dim(legs.l1)+log2dim(legs.l2))  # only counts external legs
end

function greedy_loss(::MinSpaceDiff, incidence_list, log2_edge_sizes, vi, vj)
    log2dim(legs) = sum(l->log2_edge_sizes[l], legs, init=0)
    legs = analyze_contraction(incidence_list, vi, vj)
    D1,D2,D12,D01,D02,D012 = log2dim.(getfield.(Ref(legs), 1:6))
    exp2(D01+D02+D012) - exp2(D1+D01+D12) - exp2(D2+D02+D12)  # out - in
end

function space_complexity(incidence_list, log2_sizes)
    sc = 0.0
    for v in vertices(incidence_list)
        for e in edges(incidence_list, v)
            sc += log2_sizes[e]
        end
    end
    return sc
end