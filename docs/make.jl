using Documenter, FactorGraph

makedocs(;
    modules=[FactorGraph],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/GiggleLiu/FactorGraph.jl/blob/{commit}{path}#L{line}",
    sitename="FactorGraph.jl",
    authors="JinGuo Liu",
    assets=String[],
)

deploydocs(;
    repo="github.com/GiggleLiu/FactorGraph.jl",
)
