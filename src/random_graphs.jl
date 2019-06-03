export demo_fg, demo_fg2, random_simple_fg, random_fg
using StatsBase

function demo_fg()
    orders = [1, 3, 3, 3, 2, 4, 2]
    tensors = [randn(repeat([2], o)...) for o in orders]
    lm = zeros(Int, 7, 9)
    lm[1,1] = 1

    lm[2,1] = 1
    lm[2,2] = 2
    lm[2,3] = 3

    lm[3,2] = 1
    lm[3,4] = 2
    lm[3,5] = 3

    lm[4,3] = 1
    lm[4,6] = 2
    lm[4,7] = 3

    lm[5,4] = 1
    lm[5,8] = 2

    lm[6,5] = 1
    lm[6,6] = 2
    lm[6,8] = 3
    lm[6,9] = 4

    lm[7,7] = 1
    lm[7,9] = 2
    FGraph(tensors, lm)
end

# non-simple graph
function demo_fg2()
    orders = [2, 4, 4, 2, 2, 4, 2]
    tensors = [randn(repeat([2], o)...) for o in orders]
    lm = zeros(Int, 7, 9)
    lm[1,1] = 1

    lm[2,1] = 1
    lm[2,2] = 2
    lm[2,3] = 3
    lm[2,6] = 4

    lm[3,2] = 1
    lm[3,4] = 2
    lm[3,5] = 3
    lm[3,6] = 4

    lm[4,3] = 1
    lm[4,7] = 2

    lm[5,4] = 1
    lm[5,8] = 2

    lm[6,5] = 1
    lm[6,8] = 3
    lm[6,7] = 4

    lm[7,7] = 2
    lm[7,9] = 1
    FGraph(tensors, lm)
end

function random_simple_fg(::Type{T}, tensor_sizes::Vector{<:Tuple}, nbond::Int; bias_factor::Real=-1) where T
    n = length(tensor_sizes)
    nlegs = length.(tensor_sizes) |> sum
    nlegs < 2nbond && throw(ArgumentError("Number of bond $nbond too large."))
    lm = zeros(Int, n, nbond)
    # randomly connect tensors
    # assign tensor legs
    #lm[rand(n, nbond)<bond_density] .= 1

    # combine tensors with graph
    ts = [randn(T, size...) for size in tensor_sizes]
    fg = FGraph(ts, lm)

    for i=1:100
        try
            for ibond = 1:nbond
                legs = fg |> dangling_legs
                it, jt = sample(1:n, Weights(exp.(-bias_factor*length.(legs))), 2, replace=false)
                ileg = sample(dangling_legs(fg, it))
                jleg = sample(dangling_legs(fg, jt))
                legmap(fg)[it, ibond] = ileg
                legmap(fg)[jt, ibond] = jleg
            end
            break
        catch
            legmap(fg) .= 0
        end
        if i==100
            throw(ArgumentError("Can not construct this random simple graph with 100 tries, try decrease `bias_factor`."))
        end
    end
    fg
end

function random_fg(::Type{T}, nv::Int, ne::Int, link_density::Real=2/nv, bond_dimension::Int=2) where T
    lm = zeros(Int, nv, ne)
    # randomly connect tensors
    ks = zeros(Int, nv)
    for iv = 1:nv
        for ie = 1:ne
            if rand() < link_density
                ks[iv] += 1
                @inbounds lm[iv, ie] = ks[iv]
            end
        end
    end

    # combine tensors with graph
    ts = [randn(T, fill(bond_dimension, ndim)...) |> asarray for ndim in ks]
    FGraph(ts, lm)
end

random_fg(nv::Int, ne::Int, args...) = random_fg(Float64, nv::Int, ne::Int, args...)
random_simple_fg(tensor_sizes::Vector{<:Tuple}, nbond::Int; kwargs...) = random_simple_fg(Float64, tensor_sizes, nbond; kwargs...)
