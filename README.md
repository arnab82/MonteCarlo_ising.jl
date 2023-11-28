# Monte Carlo Simulation for 2D Ising Model
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/arnab82/MonteCarlo_ising.jl/.github/workflows/blank.yml//badge.svg)

## Overview

This repository contains a Julia implementation of a Monte Carlo simulation for the 2D Ising model. The Ising model is a mathematical model in statistical mechanics that describes the magnetic properties of a two-dimensional lattice of spins.

## Features

- Monte Carlo simulation of a 2D Ising model.
- Visualization of the spin configurations.
- Calculation of thermodynamic properties such as magnetization, energy, and heat capacity.

## Installation
```bash
git clone https://github.com/arnab82/MonteCarlo_ising.jl.git
cd MonteCarlo_ising.jl
julia --project=./ -tauto
using Pkg
Pkg.instantiate()
Pkg.add("Plots")
Pkg.precompile()
using MonteCarlo_ising
Pkg.test()
```
## Copyright

Copyright (c) 2023, Arnab Bachhar
