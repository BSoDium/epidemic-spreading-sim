"""
Create a new state for a contact network and infect certain ammount of nodes.

**PARAMS** :
  - nv `Int` : *number of vertices*
  - ni `Int` : *number of infected*
**RETURNS** :
  - state `Array{Int32,1}` : *The newly created state of the contact network*
  
"""
function init_State(nv, ni)
  # Initialize state array
  state = [0 for i in 0:nv]
  
  # Infect ni nodes (randomly chosen)
  UninfectedIndices = [i for i in 1:nv]
  for i = 1:min(ni,nv)
    randi = rand(UninfectedIndices)
    UninfectedIndices = setdiff(UninfectedIndices, [randi])
    state[randi] = 1
  end

  return state
end

"""
Draw the contact network and the state of its vertices.

**PARAMS** :
  - net `LightGraph` : *graph representing the contact network*
  - state `Array{Int32,1}` : *disease status of each vertex*
  - filename `String` : *name of the file to save the image*
  - susceptibleColor `String` : *color of the susceptible vertices (optional)*
  - infectedColor `String` : *color of the infected vertices (optional)*
  - recoveredColor `String` : *color of the recovered vertices (optional)*
  - cleanColor `String` : *color of the clean vertices (optional)*
  
"""
function draw_Net_Graph(
    net, 
    state, 
    filename, 
    susceptibleColor = colorant"lightseagreen", 
    infectedColor = colorant"orange", 
    recoveredColor = colorant"lightgreen", 
    alertColor = colorant"purple"
  )
  nodecolor = [susceptibleColor for i in 0:nv(net)]
  for i in 1:nv(net)
    if state[i] == 1
      nodecolor[i] = infectedColor
    elseif state[i] == 2
      nodecolor[i] = recoveredColor
    elseif state[i] == 3
      nodecolor[i] = alertColor
    end
  end
  draw(PNG(string("out/", filename), 100cm, 100cm), gplot(net, nodefillc=nodecolor))
  # gplot(net, nodefillc=nodecolor)
end