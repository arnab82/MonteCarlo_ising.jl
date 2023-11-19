using MonteCarlo_ising
using JLD2

lattice_n =MonteCarlo_ising.IsingModel(50, 1.0, true)
display(lattice_n)
enrg = MonteCarlo_ising.energy_manual(lattice_n)
display(enrg)
@save "lattice_n.jld2" lattice_n enrg