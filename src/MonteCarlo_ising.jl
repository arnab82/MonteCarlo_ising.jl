module MonteCarlo_ising


using LinearAlgebra
using Random
using Printf
using Statistics
using Plots
include("Ising_model.jl")
include("metropolis.jl")
include("properties.jl")

end # module
