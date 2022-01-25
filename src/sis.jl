using Pkg
using LightGraphs
using GraphPlot
using Colors
using CairoMakie
using StatsBase
using Plots
using JLD2
using Compose
using Printf

include("lib.jl")

"""
Take a contact network at a certain state and apply t time steps
of a SIS model.

**PARAMS** :
  - net `LightGraph` : *graph representing the contact network*
  - state `Array{Int32,1}` : *disease status of each vertex*
  - beta `Float64` : *infection rate*
  - alpha `Float64` : *curing rate*
  - t `Int32` : *number of time step*

**RETURNS** :
  - state `Array{Int32,1}` : *The new state of the contact network after t time steps.*
  - infection_percentage `Array{Float64,1}` : *The average percentage of infected nodes at each time step.*
"""
function Core_SIS(net,state,beta,alpha,t)
  
  # This array will contain the average infection percentage through each time
  # step.
  infection_percentage = zeros(t)

  for i in 1:t
    infected_cells = LightGraphs.findall(==(1), state)
    for cell in infected_cells
      neighbors = LightGraphs.neighbors(net, cell)
      for neighbor in neighbors
        # infect neighbors
        if state[neighbor] == 0
          if rand() < beta
            state[neighbor] = 1
          end
        end

        # cure itself
        if rand() < alpha
          state[cell] = 0
        end
      end
    end
    infection_percentage[i] = sum(state)/nv(net)
  end

  return state, infection_percentage
end

"""
Take a contact network at a certain state and apply t time steps
of a SIS model.

**PARAMS** :
  - net `LightGraph` : *graph representing the contact network*
  - state `Array{Int32,1}` : *disease status of each vertex*
  - beta `Float64` : *infection rate*
  - alpha `Float64` : *curing rate*
  - t `Int32` : *number of time step*

**RETURNS** :
  - state `Array{Int32,1}` : *The new state of the contact network after t time steps.*
  
"""
function SIS(net, state, beta, alpha, t)
  return Core_SIS(net, state, beta, alpha, t)[1]
end

"""
Take a contact network, different diseases (defined by
different parameters alpha and beta), a number of initial
infected people and process nbsimu simulations of SIS over
t time steps. You will provide the prediction of the
percentage of infected at each time t as well as the
spreading rate of each disease.

**PARAMS** :
  - net `LightGraph` : *graph representing the contact network*
  - nbinf `Int32` : *number of infected at the start of each*
        simulation
  - betas `Array{Float64,1}` : *array of infection rate on edges*
  - alphas `Array{Float64,1}` : *array of curing rate on vertices*
  - t `Int32` : *number of time steps*
  - nbsimu `Int32` : *number of simulations*

**RETURNS** :
  - `Array{Float64,2}` : *the prediction of the percentage of 
        infected at each time step and for each disease. The 
        first dimension contains the time steps and the second
        contains the diseases*
  - `Array{Float64,1}` : *effective spreading rate for each 
        disease*

"""
function Simulation_SIS(net,nbinf,betas,alphas,t,nbsimu)
  # initialize lock
  lk = Threads.SpinLock()

  nbdis = size(alphas)[1] 
  avg_infection_percentage = zeros(t, nbdis)
  effective_spreading_rate = zeros(nbdis)

  Threads.@threads for i in 1:nbdis
    alpha = alphas[i]
    beta = betas[i]
    effective_spreading_rate[i] = beta / alpha
    for j in 1:nbsimu
      # Initialize the state
      state = init_State(nv(net), nbinf)
      # Simulate the disease
      state, infection_percentage = Core_SIS(net, state, beta, alpha, t)

      lock(lk)
      try
        # Compute the percentage of infected
        avg_infection_percentage[: , i] += infection_percentage
      finally
        unlock(lk)
      end
    end
  end
  return avg_infection_percentage / nbsimu , effective_spreading_rate
end

if abspath(PROGRAM_FILE) == @__FILE__
  # just a small test
  karat7 = smallgraph(:karate)

  betas=[0.05,0.1,0.01,0.4,0.04,0.05,0.005]
  alphas=[0.05,0.1,0.01,0.1,0.01,0.1,0.01]

  predictions, taus = Simulation_SIS(karat7,2,betas,alphas,10,7)

  Plots.plot(predictions, label=taus',xlabel="Time",ylabel="Avg % of infected")
end