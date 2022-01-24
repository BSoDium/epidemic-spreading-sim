"""
Take a contact network at a certain state and apply t time steps
of an SAIR model.

**PARAMS**
  - net `LightGraph`: *graph representing the contact network*
  - state `Array{Int32,1}`: *disease status of each vertex*
  - beta0 `Float64`: *infection rate when not alert*
  - beta1 `Float64`: *infection rate when alert*
  - alpha `Float64`: *curing rate*
  - kappa `Float64`: *alerting rate*
  - t `Int32`: *number of time steps*

**RETURNS**
  - `Array{Int32,1}`: *The new state of the contact network after t time steps.*
"""      
function SAIR(net,state,beta0,beta1,alpha,kappa,t)
  # This array will let us know if we have already checked if the node should
  # cured at this step or not.
  hasCellBeenChecked = [0 for i in 0:nv(net)]

  # Initialize the new state
  outputstate = state
  for i in 1:t
      for j in LightGraphs.edges(net)
          from = j.src
          to = j.dst

          # check if the nodes can infect each other
          # not alert
          if state[from] == 1 && state[to] == 0
              if rand() < beta0
                  outputstate[to] = 1
              end
          elseif state[to] == 1 && state[from] == 0
              if rand() < beta0
                  outputstate[from] = 1
              end
          end

          # alert
          if state[from] == 1 && state[to] == 3
              if rand() < beta1
                  outputstate[to] = 1
              end
          elseif state[to] == 1 && state[from] == 3
              if rand() < beta1
                  outputstate[from] = 1
              end
          end

          # check if the infected nodes can alert the susceptible nodes
          if state[from] == 1 && state[to] == 0
              if rand() < kappa
                  outputstate[to] = 3
              end
          elseif state[to] == 1 && state[from] == 0
              if rand() < kappa
                  outputstate[from] = 3
              end
          end


          # check both cells can cure themselves, if they have not been checked yet
          if state[from] == 1 && hasCellBeenChecked[from] == 0
              if rand() < alpha
                  outputstate[from] = 2
              end
              hasCellBeenChecked[from] = 1
          elseif state[to] == 1 && hasCellBeenChecked[to] == 0
              if rand() < alpha
                  outputstate[to] = 2
              end
              hasCellBeenChecked[to] = 1
          end
      end
  end

  return outputstate
end

"""
Take a contact network, different diseases (defined by 
different parameters alpha and beta), a number of initial
infected people and process nbsimu simulations of SAIR over
t time steps. You will provide the prediction of the 
percentage of infected at each time t as well as the 
spreading rate of each disease.

**PARAMS**
  - net `LightGraph`: *graph representing the contact network*
  - nbinf `Int32`: *number of infected at the start of each 
      simulation*
  - betas0 `Array{Float64,1}`: *array of infection rate when not alert on edges*
  - betas1 `Array{Float64,1}`: *array of infection rate when alert on edges*
  - alphas `Array{Float64,1}`: *array of curing rate on vertices*
  - kappas `Array{Float64,1}`: *array of alerting rate on edges*
  - t `Int32`: *number of time step*
  - nbsimu `Int32`: *number of simulations*

**RETURNS**
  - `Array{Float64,3}`: *the prediction of the percentage of 
      infected, the percentage of susceptible and the 
      percentage of recovered at each time step and for each 
      disease. The first dimension contains the time steps,
      the second contains the diseases, and the third the status
      (Infected: `[:,:,1]`, Recovered: `[:,:,2]`, Susceptible: `[:,:,3]`)*
  - `Array{Float64,1}`: *effective spreading rate for each 
      disease*
      
"""
function Simulation_SAIR(net,nbinf,betas0,betas1,alphas,kappas,t,nbsimu)
  
# initialize lock
lk = Threads.SpinLock()

nbdis = size(alphas)[1]
avg_infected_percentage = zeros(t, nbdis, 3)
effective_spreading_rate = zeros(nbdis)

Threads.@threads for i in 1:nbdis
alpha = alphas[i]
beta0 = betas0[i]
beta1 = betas1[i]
kappa = kappas[i]

effective_spreading_rate[i] = beta0 / alpha
for j in 1:nbsimu
    # Initialize the state
    state = init_State(nv(net), nbinf)
    for step in 1:t
    # Simulate the disease
    state = SAIR(net, state, beta0, beta1, alpha, kappa, step)

    lock(lk)
    try
        # Compute the percentage of infected
        for s in state
        if s == 1
            avg_infected_percentage[step, i, 1] += 1
        elseif s == 2
            avg_infected_percentage[step, i, 2] += 1
        elseif s == 3
            avg_infected_percentage[step, i, 3] += 1
        end
        end
        avg_infected_percentage[step, i, 1] /= nv(net)
        avg_infected_percentage[step, i, 2] /= nv(net)
        avg_infected_percentage[step, i, 3] /= nv(net)
    finally
        unlock(lk)
    end
    end
end
end
return avg_infected_percentage / nbsimu , effective_spreading_rate
end
end