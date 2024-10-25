using QuantumCollocationCore
using Documenter

DocMeta.setdocmeta!(QuantumCollocationCore, :DocTestSetup, :(using QuantumCollocationCore); recursive=true)

makedocs(;
    modules=[QuantumCollocationCore],
    authors="Aaron Trowbridge <aaron.j.trowbridge@gmail.com> and contributors",
    sitename="QuantumCollocationCore.jl",
    format=Documenter.HTML(;
        canonical="https://aarontrowbridge.github.io/QuantumCollocationCore.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aarontrowbridge/QuantumCollocationCore.jl",
    devbranch="main",
)
