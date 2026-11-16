# Set up the Julia reference environment used to generate test fixtures and
# validate SNaQ++ against SNaQ.jl. Expects a checkout of SNaQ.jl (the branch
# being ported) at ../../SNaQ.jl relative to this file.
using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(path=joinpath(@__DIR__, "..", "..", "SNaQ.jl"))
Pkg.instantiate()
Pkg.precompile()
