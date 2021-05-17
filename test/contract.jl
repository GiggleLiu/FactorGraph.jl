using FactorGraph
using OMEinsum
using Zygote
using Test

function bpcheck(f, args...; η = 1e-5, verbose = false)
    g = gradient(f, args...)
    all(x->x===(nothing), g) && error()
    dy_ref = 0
    for x in g
        x === nothing && continue
        x isa Tuple && (dy_ref += η * mapreduce(y -> y === nothing ? 0 : sum(abs2,y), +, x))
        x isa AbstractArray && (dy_ref += η * sum(abs2,x))
    end
    dy = f(args...) - f([gi === nothing ? arg : arg .- η .* gi for (arg, gi) in zip(args,g)]...)

    verbose && @show dy
    verbose && @show dy_ref

    isapprox(dy, dy_ref, rtol=1e-2, atol=1e-8)
end

#Base.conj(arr::AbstractArray{<:Any, 0}) = conj!(copy(arr))

@testset "einsum bp" begin
    a = randn(ComplexF64, 3,3)
    f2(a) = ein"ij,ik,il->jkl"(a, conj(a), a) |> sum |> real

    @test gradient(f2, a)[1] !== nothing
    @test bpcheck(f2, a)
    @test bpcheck(a->treecontract(((1,2),3), ((1,2), (2,3), (3,1)), (a, a, a), ())[] |> real, a)
end

@testset "optcontract" begin
    a = randn(ComplexF64, 3,3)
    @test optcontract(((1,2), (2,3)), (a, a), (1,3)) ≈ a*a
    @test bpcheck(a->ein"ij,ji->"(a, a) |> sum |> real, a)
    @test bpcheck(a->optcontract(((1,2), (2,3)), (a, a), (1,3)) |> sum |> real, a)
    @test bpcheck(a->optcontract(((1,2), (2,3), (3,4)), (a,a,a), (4,1)) |> sum |> real, a)
end