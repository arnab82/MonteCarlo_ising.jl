using MonteCarlo_ising
using JLD2
using Test
function run()
    @load "lattice_n.jld2" 
    @test MonteCarlo_ising.energy_manual(lattice_n) == enrg
end
run()
