# QuantumCollocationCore

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kestrelquantum.github.io/QuantumCollocationCore.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kestrelquantum.github.io/QuantumCollocationCore.jl/dev/)
[![Build Status](https://github.com/kestrelquantum/QuantumCollocationCore.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kestrelquantum/QuantumCollocationCore.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kestrelquantum/QuantumCollocationCore.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kestrelquantum/QuantumCollocationCore.jl)


This package provides a core library for quantum collocation methods. It is designed to be used in conjunction with the [QuantumCollocation.jl](https://github.com/kestrelquantum/QuantumCollocation.jl) package, which provides a high-level interface for solving quantum optimal control problems using direct collocation.

The underlying nonlinear solver is [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl), which is a Julia interface to the [Ipopt](https://coin-or.github.io/Ipopt/) solver. 

## Quickstart developers guide

__Install Julia__ [Juliaup](https://github.com/JuliaLang/juliaup) is an installer and version manager. This is one useful way to manage Julia versions and keep up with the latest changes. After installing, run `julia` to obtain the Julia _REPL_.

__Julia environments__
[(Documentation)](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project) Your project's environment is stored in _Project.toml_. You can interactively add packages to an environment by using the Julia command line _REPL_ and _package manager_.  Start Julia in the project folder. Type `]` to enter the package manager. Type `activate .` to activate or create an environment specified by _Project.toml_ located in the current folder. Separately, you generate a manifest (solving the versions to create a valid environment) by running `instantiate`; instantiate will check that the environment is correct after you add all the packages you want.

__Adding packages__
[(Documentation)](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-packages) The initial cell for a Piccolo notebook might look something like the following:
```Julia
# Standard packages
using LinearAlgebra
using CairoMakie

# Piccolo packages
using QuantumCollocationCore
using NamedTrajectories
using TrajectoryIndexingUtils
```

First, let's install some standard packages (these are like Numpy and Matplotlib). Open the package manager in the current environment (type `julia`, `]`, and `activate .`), type `add LinearAlgebra` to install and precompile _LinearAlgebra_. Same with `CairoMakie`. 

Second, let's install _Piccolo_. There are three packages (_QuantumCollocation_, _NamedTrajetories_, _TrajectoryIndexingUtils_) inside [Piccolo](https://docs.juliahub.com/General/Piccolo/stable/). We could do `add Piccolo` to get the three as a bundle from the Julia repository. Instead of individually calling `using ...` for each, this approach only requires `using Piccolo` at the start of a file or notebook.

As a developer, we want to use the git repositories directly from [the Kestrel Quantum Github page](https://github.com/kestrelquantum). Clone, then add the local packages to the Project file with e.g. `dev ../relative/path/to/repo/QuantumCollocation`. This command installs the development version of _QuantumCollocation_ pointing to the local Github code instead of the package repository. You can repeat this for the others, also.

__Developing__
[Revise.jl](https://timholy.github.io/Revise.jl/stable/) will let you edit source code, update packages, and reload the changes in a notebook---automatically! This is a great tool for development. `add Revise` from the REPL and then include it before any packages you intend to edit:
```Julia
using Revise
using QuantumCollocationCore
```

### Documentation

Documentation is built using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) and uses [Literate.jl](https://github.com/fredrikekre/Literate.jl) to generate markdown files from scripts stored in [docs/literate](docs/literate). To build the documentation locally, start julia with the docs environment:

```bash
julia --project=docs
```

Then run the following commands in the Julia REPL:

```
] instantiate
```

Then (for ease of development) load the following packages:

```julia
using Revise, LiveServer, QuantumCollocation
```

To live-serve the docs, run
```julia
servedocs(literate_dir="docs/literate", skip_dir="docs/src/generated")
```

Changes made to files in the docs directory should be automatically reflected in the live server. To reflect changes in the source code (e.g. doc strings), since we are using Revise, simply kill the live server running in the REPL (with, e.g., Ctrl-C) and restart it with the above command. 

### Tips for Visual Studio Code
__Julia extension__ You can run Julia notebooks and much more with [the Julia extension](https://code.visualstudio.com/docs/languages/julia). Upon opening your project folder in VS code and attempting to run an `.ipynb`, you will see that VS Code finds the interpreters managed by juliaup and defaults to using the environment based on the _Project.toml_ in the project directory.

__Fonts__ VS Code will not display all characters allowed by Julia. You can change the editor font family in the settings to `'JuliaMono'` to get full support. If you don't want to mix and mash, you can create a new VS Code settings profile for working in Julia at _File>Preferences>Profile_.

__Tests__ Tests should automatically populate in VS Code when working with a Piccolo package. For example, just by adding the `QuantumCollocation.jl` folder to your workspace, you should see tests appear if you click on the _Testing_ sidebar icon. If you run one of these tests, a new Julia kernel is spawned for the test. You can find the kernel if you click on the _Julia_ sidebar icon (after installing the Julia extensions). Sometimes, for the tests to recognize new changes, you may need to manually kill this kernel to see your changes reflected.

