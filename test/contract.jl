using FactorGraph
using OMEinsum: bpcheck, einsum
using Zygote
using Test

Base.conj(arr::AbstractArray{<:Any, 0}) = conj!(copy(arr))

@testset "einsum bp" begin
    a = randn(ComplexF64, 3,3)
    f2(a) = einsum(((1,2), (1,3), (1,4)), (a, conj(a), a), (2,3,4)) |> sum |> real

    @test gradient(f2, a)[1] !== nothing
    @test bpcheck(f2, a)
    @test bpcheck(a->treecontract(((1,2),3), ((1,2), (2,3), (3,1)), (a, a, a), ())[] |> real, a)
end

@testset "optcontract" begin
    a = randn(ComplexF64, 3,3)
    @test optcontract(((1,2), (2,3)), (a, a), (1,3)) ≈ a*a
    @test bpcheck(a->einsum(((1,2), (2,1)), (a, a), ()) |> sum |> real, a)
    @test bpcheck(a->optcontract(((1,2), (2,3)), (a, a), (1,3)) |> sum |> real, a)
    @test bpcheck(a->optcontract(((1,2), (2,3), (3,4)), (a,a,a), (4,1)) |> sum |> real, a)
end
