"""
Take a contact network at a certain state and apply t time steps
  of an SIR model.
  
  PARAMS
    - net `LightGraph`: *graph representing the contact network*
    - state `Array{Int32,1}`: *disease status of each vertex*
    - beta `Float64`: *infection rate*
    - alpha `Float64`: *curing rate*
    - t `Int32`: *number of time step*
  
  RETURNS
    - outputstate `Array{Int32,1}`: *The new state of the contact network after t time steps.*
    - infection_percentage `Array{Float64,1}` : *The average percentage of infected nodes at each time step.*
  """
function Core_SIR(net,state,beta,alpha,t)

  stats = zeros(t,3);

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
          state[cell] = 2
        end
      end
    end
    step_stats = [
      LightGraphs.count(==(1), state),
      LightGraphs.count(==(2), state),
      LightGraphs.count(==(0), state)
    ]
    stats[i, :] = step_stats
  end

  return state, stats
end

"""
Take a contact network at a certain state and apply t time steps
  of an SIR model.
  
  PARAMS
    - net `LightGraph`: *graph representing the contact network*
    - state `Array{Int32,1}`: *disease status of each vertex*
    - beta `Float64`: *infection rate*
    - alpha `Float64`: *curing rate*
    - t `Int32`: *number of time step*
  
  RETURNS
      - `Array{Int32,1}`: *The new state of the contact network after t time steps.*
  """
function SIR(net,state,beta,alpha,t)
  return Core_SIR(net,state,beta,alpha,t)[1]
end


"""
Take a contact network, different diseases (defined by 
different parameters alpha and beta), a number of initial
infected people and process nbsimu simulations of SIR over
t time steps. You will provide the prediction of the 
percentage of infected at each time t as well as the 
spreading rate of each disease.

**PARAMS**
  - net `LightGraph` : *graph representing the contact network*
  - nbinf `Int32` : *number of infected at the start of each 
      simulation*
  - betas `Array{Float64,1}` : *array of infection rate on edges*
  - alphas `Array{Float64,1}` : *array of curing rate on vertices*
  - t `Int32` : *number of time steps*
  - nbsimu `Int32` : *number of simulations*

RETURNS
  - `Array{Float64,3}` : *the prediction of the percentage of 
      infected, the percentage of susceptible and the 
      percentage of recovered at each time step and for each 
      disease. The first dimension contains the time steps,
      the second contains the diseases, and the third one of the
      following statuses:
      - Infected : `[:, :, 1]`
      - Recovered : `[:, :, 2]`
      - Susceptible : `[:, :, 3]`
  - `Array{Float64,1}` : *effective spreading rate for each 
      disease*

"""
function Simulation_SIR(net,nbinf,betas,alphas,t,nbsimu)
  # initialize lock
  lk = Threads.SpinLock()

  nbdis = size(alphas)[1]
  avg_infection_stats = zeros(t, nbdis, 3)
  effective_spreading_rate = zeros(nbdis)

  for i in 1:nbdis
    alpha = alphas[i]
    beta = betas[i]
    effective_spreading_rate[i] = beta / alpha
    for j in 1:nbsimu
      # Initialize the state
      state = init_State(nv(net), nbinf)
        # Simulate the disease
        state, stats = Core_SIR(net, state, beta, alpha, t)

        lock(lk)
        try
          avg_infection_stats[:, i, 1] += stats[:, 1]
        finally
          unlock(lk)
        end
      end
    end
  return avg_infection_stats / nbsimu , effective_spreading_rate
end