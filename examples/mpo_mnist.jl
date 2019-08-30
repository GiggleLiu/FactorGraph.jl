using Zygote, OMEinsum
using Zygote:relu
using Flux, Flux.Data.MNIST
using Statistics
using Flux.Optimise

imgs = MNIST.images()

# Stack images into one large batch
X = hcat(float.(reshape.(imgs, :))...) |> gpu

labels = MNIST.labels()
# One-hot-encode the labels
Y = Flux.onehotbatch(labels, 0:9)

struct MPO{T}
  tensors::Vector{Array{T,4}}
end

Base.Matrix(mps::MPO) = reshape(ein"aijb,bklc,cmnd,dope,eqrf->ikmoqjlnpr"(mps.tensors...), 32, 32)

function loss(x, W1, W2, W3, y)
    y1 = relu.(ein"ij,ki->kj"(x, W1))
    y2 = relu.(ein"ij,ki->kj"(y1, Matrix(W2)))
    q = softmax(ein"ij,ki->kj"(y2, W3))
    #mean(Flux.crossentropy(y, q)))
    mean(abs2.(y .- q))
end

_kldivergence(p, q) = p*log(p/q)

function pick_batch(X, Y, nbatch::Int=100)
    picked = rand(1:60000, nbatch)
    X[:,picked], Y[:,picked]
end

function train(X, Y; optimizer, niter)
    x, y = pick_batch(X, Y)
    W1 = randn(32, 784)
    W3 = randn(10, 32)
    mpo = MPO([randn(1,2,2,4), randn(4,2,2,4), randn(4,2,2,4), randn(4,2,2,4), randn(4,2,2,1)])

    for i=1:niter
        gx, gW1, gmpo, gW3, _ = Zygote.gradient(loss, x, W1, mpo, W3, y)
        for (data, gdata) in zip((x, W1, W3), (gx, gW1, gW3))
            Optimise.update!(optimizer, data, gdata)
        end
        for (data, gdata) in zip(mpo.tensors, gmpo.tensors)
            Optimise.update!(optimizer, data, gdata)
        end
        @show loss(x, W1, mpo, W3, y)
    end
end

optimizer = Flux.Optimise.ADAM(0.01)
train(X, Y; optimizer=optimizer, niter=100)
