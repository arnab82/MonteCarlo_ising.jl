
function magnetization(spin_array::IsingModel, skipped_points::Int64, total_points::Int64, T::Float64, k::Float64=1.0)
    J = spin_array.J
    bj = J / (k * T)
    lattice = spin_array.lattice
    N = spin_array.N

    spin, energies, spins = metropolis(spin_array, total_points, bj, energy_manual(spin_array))
    ms_j = abs(mean(spin[skipped_points:end]))
    ms_squared_j = mean(spin[skipped_points:end] .^ 2)
    return ms_j, ms_squared_j
end

function magnetic_susceptibility(spin_array::IsingModel, skipped_points::Int64, total_points::Int64, T)
    Ms, Ms_squared = magnetization(spin_array, skipped_points, total_points, T)
    N = spin_array.N
    k = 1.0
    bj = spin_array.J / (k * T)
    T_c = 2.269 * spin_array.J / k
    if T > T_c
        chi = (Ms_squared) * (N^2) / (k * T)
    elseif T < T_c
        chi = (Ms_squared - Ms^2) * (N^2) / (k * T)
    end

    return chi
end

function standard_deviation_spin(spin_array::IsingModel, skipped_points::Int64, total_points::Int64, T) 
    Ms, Ms_squared = magnetization(spin_array, skipped_points, total_points, T)
    return sqrt(Ms_squared - Ms^2)
end
"""
# Heat Capacity
Heat capacity, C_v=frac{sigma^2_E}{T^2}
=(<E^2>-<E>^2).beta^2k^2
=sigma_{E/J}^2.(beta J)^2k^2
"""
function specific_heat(spin_array::IsingModel, skipped_points::Int64, total_points::Int64, T)
    N = spin_array.N
    k = 1.0
    bj = spin_array.J / (k * T)
    E_mean = mean(energies[skipped_points:end])
    E_squared_mean = mean(energies[skipped_points:end] .^ 2)
    return (E_squared_mean - E_mean^2) * (k * bj^2)

end
# temperatures = [1.6, 2.0, 2.16, 2.34, 2.43, 2.76, 3.35]
# skipped_points = 10000
# total_points = 200000
# plot()
# ms_values = Float64[]
# chi_values = Float64[]
# for temperature in temperatures
#     spin_array = IsingModel(50)
#     chi = magnetic_susceptibility(spin_array, skipped_points, total_points, temperature)
#     ms, ms_squared = magnetization(spin_array, skipped_points, total_points, temperature)
#     push!(ms_values, ms)
#     push!(chi_values, chi)
# end

# # Plot results for each lattice size
# plot!(temperatures, chi_values, label="L = 50")
