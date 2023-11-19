
using Plots
function periodic_boundary_conditions(x, y, N::Int64)
    """
    parameters:
    x: x coordinate of the spin
    y: y coordinate of the spin
    N: size of the lattice
    _________________________
    returns:
    x_p1: x coordinate of the spin to the right
    x_m1: x coordinate of the spin to the left
    y_p1: y coordinate of the spin to the top
    y_m1: y coordinate of the spin to the bottom
    """
    x_p1 = mod1(x + 1, N)
    x_m1 = mod1(x - 1, N)
    y_p1 = mod1(y + 1, N)
    y_m1 = mod1(y - 1, N)
    return x_p1, x_m1, y_p1, y_m1
end
function periodic_boundary_conditions_manual(x, y, N::Int64)
    """
    parameters:
    x: x coordinate of the spin
    y: y coordinate of the spin
    N: size of the lattice
    _________________________
    returns:
    x_p1: x coordinate of the spin to the right
    x_m1: x coordinate of the spin to the left
    y_p1: y coordinate of the spin to the top
    y_m1: y coordinate of the spin to the bottom
    """
    x_p1 = if x == N
        1
    else
        x + 1
    end

    x_m1 = if x == 1
        N
    else
        x - 1
    end

    y_p1 = if y == N
        1
    else
        y + 1
    end

    y_m1 = if y == 1
        N
    else
        y - 1
    end

    return x_p1, x_m1, y_p1, y_m1
end

function calculate_probability(delta_E::Float64, Bj::Float64)
    """
    parameters
    delta_E: energy difference between the initial and final state
    Bj: Boltzmann constant times temperature
    _________________________
    returns:
    probability of changing state
    """
    if delta_E <= 0
        return 1.0
    else
        return exp(-Bj * delta_E)
    end
end
function metropolis(spin_array::IsingModel, times::Int64, Bj::Float64, energy::Float64)
    """
    parameters:
    spin_array: IsingModel object
    times: number of times to run the simulation
    Bj: Boltzmann constant times temperature
    energy: initial energy of the system
    _________________________
    returns:
    net_spin: 1D array of net magnetization array 
    net_energy: 1D array of net energy 
    
    """
    N=spin_array.N
    J=spin_array.J
    spin_array = copy(spin_array.lattice)
    net_spin = zeros(times - 1)
    net_energy = zeros(times - 1)
    
    for t in 1:times - 1
        # choosing random spin and flipping it
        x = rand(1:N)
        y = rand(1:N)
        spin_i = spin_array[x, y]  # initial spin
        spin_f = -spin_i  # proposed final spin

        # calculating energy difference
        E_i = 0
        E_f = 0

        # periodic boundary conditions

        # x_p1, x_m1, y_p1, y_m1 = periodic_boundary_conditions(x, y, N)
        x_p1, x_m1, y_p1, y_m1 = periodic_boundary_conditions_manual(x, y, N)

        #calculating energy difference for spin configuration  i and f
        E_i += J*-spin_i * spin_array[x_m1, y]
        E_f += J*-spin_f * spin_array[x_m1, y]

        E_i += J*-spin_i * spin_array[x_p1, y]
        E_f += J*-spin_f * spin_array[x_p1, y]

        E_i += J*-spin_i * spin_array[x, y_m1]
        E_f += J*-spin_f * spin_array[x, y_m1]

        E_i += J*-spin_i * spin_array[x, y_p1]
        E_f += J*-spin_f * spin_array[x, y_p1]

        delta_E = E_f - E_i

        # changing state with probability

        prob = calculate_probability(delta_E, Bj)

        if rand() < prob
            spin_array[x, y] = spin_f
            energy += delta_E
        end

        net_spin[t] = sum(spin_array)
        net_energy[t] = energy

    end

    return net_spin, net_energy
end
# #For 75% negative spin configuration
# spins,energies=metropolis(lattice_n,100000,0.86,energy_manual(lattice_n))
# display(spins)
# plot_spins = plot(spins / N^2, xlabel="Time steps", ylabel="Average spin, \$\\bar{m}\$", grid=true,label="75% Negative initial spin",dpi=800)
# savefig("./metropolis_plot_s_initial_negative_spin.png")
# plot_energies = plot(energies, xlabel="Time steps", ylabel="Energy, \$E\$", grid=true,label="75% Negative initial spin",dpi=800)
# savefig("./metropolis_plot_e_initial_negative_spin.png")
# plot_layout = hcat(plot_spins, plot_energies)
# gr()
# plot_layout

using LinearAlgebra
using Statistics
function get_spin_energy(lattice::IsingModel, Bjs::AbstractArray,skipped_points::Int64,total_points::Int64)
    """
    parameters:
    lattice: IsingModel object
    Bjs: array of Boltzmann constant times temperature
    _________________________
    returns:
    ms: array of average spin
    E_mean: array of average energy
    E_stds: array of standard deviation of energy
    """
    N=lattice.N
    ms = zeros(length(Bjs))
    E_mean = zeros(length(Bjs))
    E_stds = zeros(length(Bjs))
    for (j, bj) in enumerate(Bjs)
        spins, energies = metropolis(lattice, total_points, bj, energy_manual(lattice))
        ms[j] = mean(spins[skipped_points:end]) / N^2
        E_mean[j] = mean(energies[skipped_points:end])
        E_stds[j] = std(energies[skipped_points:end])
    end

    return ms, E_mean, E_stds
end

# Bjs = range(0.1, stop=1.0, length=10)
# ms_n, E_mean_n, E_stds_n = get_spin_energy((lattice_n), Bjs)
# ms_p, E_mean_p, E_stds_p = get_spin_energy((lattice_p), Bjs)

# plt = plot(size=(800, 400), legend=:bottomright, xlabel="\$(\\frac{k}{J}T\$)", ylabel="Average spin, \$\\bar{m}\$")

# plot!(1 ./ Bjs, ms_n, marker=:o, label="75% Negative initial spin")
# plot!(1 ./ Bjs, ms_p, marker=:o, label="75% Positive initial spin")

# display(plt)


# plt = plot(size=(800, 400), legend=:topright, xlabel="\$(\\frac{k}{J}T\$)", ylabel="\$(C_v/k^2\$)")

# plot!(1 ./ Bjs, E_stds_n .* Bjs, marker=:o, label="75% Negative initial spin")
# plot!(1 ./ Bjs, E_stds_p .* Bjs, marker=:o, label="75% Positive initial spin")

# display(plt)
# using Plots
# Ns = [100, 90,80, 70, 60, 50]
# J=0.8
# plt = plot(size=(800, 400), legend=:topright, xlabel="Time steps", ylabel="Average spin, \$\\bar{m}\$", grid=false)
# for N in Ns
#     lattice_n= IsingModel(N, 0.8,false)
#     spins,energies=metropolis(lattice_n,100000,0.86,energy_manual(lattice_n))
#     plot!(spins / N^2, xlabel="Time steps", ylabel="Average spin, \$\\bar{m}\$", grid=true,label="N=$N")
# end
# savefig("./metropolis_plot_negative_spin_N.png")
