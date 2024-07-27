"""
## Generating Simulation Data for Agent Based Model

To generate simulation data using the agent based model, you can use the `generate_simulation_data_agent` function defined in the `UCIWWEIHR.jl` package. This function allows you to customize various parameters for the simulation.
NOT FINISHED, STILL NEEDS WW AND RT

### Function Signature

# Arguments
- seed::Int64: Seed for random number generation. Default value is 1.
- pop_size::Int64: Size of the population. Default value is 1000.
- I_init::Int64: Initial number of infected individuals. Default value is 200.
- H_init::Int64: Initial number of hospitalized individuals. Default value is 20.
- beta::Float64: Transmission rate. Default value is 0.001.
- gamma::Float64: Rate of exposed individuals becoming infectious. Default value is 1/4.
- nu::Float64: Rate of infected individuals recovering or getting hospitalized. Default value is 1/7.
- epsilon::Float64: Rate of hospitalized individuals recovering. Default value is 1/5.
- w_init::Float64: Probability of an infected individual becoming hospitalized. Default value is 0.35.

# Returns
- df::DataFrame: A DataFrame containing the simulation data with columns `Time`, `S`, `E`, `I`, `H`, and `R`.
"""
function generate_simulation_data_agent(
    seed::Int64=1, pop_size::Int64=1000, 
    I_init::Int64=200, H_init::Int64=20, 
    beta::Float64=0.001, gamma::Float64=1/4, 
    nu::Float64=1/7, epsilon::Float64=1/5, 
    w_init::Float64=0.35
)
    
    Random.seed!(seed)

    pop_vec = repeat(["S"], pop_size)
    pop_vec[1:I_init] = repeat(["I"], I_init)
    pop_vec[I_init+1:I_init+H_init] = repeat(["H"], H_init)

    rate_vec = zeros(pop_size)
    id_vec = 1:pop_size

    t = 0.0

    rate_frame = DataFrame(id = id_vec, state = pop_vec, rate = rate_vec)

    n_suseptible = sum(rate_frame.state .== "S"); n_exposed = sum(rate_frame.state .== "E")
    n_infected = sum(rate_frame.state .== "I"); n_hospitalized = sum(rate_frame.state .== "H")
    n_recovered = sum(rate_frame.state .== "R")

    state_frame = DataFrame(T = Float64[], S = Int[], E = Int[], I = Int[], H = Int[], R = Int[])
    push!(state_frame, (t, n_suseptible, n_exposed, n_infected, n_hospitalized, n_recovered))

    rate_frame.rate .= ifelse.(rate_frame.state .== "S", beta * n_infected,
                            ifelse.(rate_frame.state .== "E", gamma,
                                ifelse.(rate_frame.state .== "I", nu,
                                    ifelse.(rate_frame.state .== "H", epsilon, 0))))

    while n_infected > 0 || n_hospitalized > 0
        #global t; 
        #global n_suseptible; global n_exposed; global n_infected; global n_hospitalized; global n_recovered
        #global state_frame = copy(state_frame); global rate_frame = copy(rate_frame)
        next_event = rand(Exponential(1/sum(rate_frame.rate)),1)

        which_id = sample(rate_frame.id, Weights(rate_frame.rate), 1)

        rate_frame.state[which_id] .= ifelse(rate_frame.state[which_id] == ["S"], "E",
                                        ifelse(rate_frame.state[which_id] == ["E"], "I",
                                            ifelse(rate_frame.state[which_id] == ["I"], ifelse.(rand() < w_init, "H", "R"), "R")))

        n_infected_old = n_infected
        n_suseptible = sum(rate_frame.state .== "S"); n_exposed = sum(rate_frame.state .== "E")
        n_infected = sum(rate_frame.state .== "I"); n_hospitalized = sum(rate_frame.state .== "H")
        n_recovered = sum(rate_frame.state .== "R")

        #new_infected = n_infected - n_infected_old
        #Rt = ?
    
        rate_frame.rate .= ifelse.(rate_frame.state .== "S", beta * n_infected,
                                ifelse.(rate_frame.state .== "E", gamma,
                                    ifelse.(rate_frame.state .== "I", nu,
                                        ifelse.(rate_frame.state .== "H", epsilon, 0))))

        t += next_event[1]
        current_state = (t, n_suseptible, n_exposed, n_infected, n_hospitalized, n_recovered)
        state_frame = push!(state_frame, current_state)
    
    end

    state_frame.Time = ceil.(Int, state_frame.T)
    state_frame.Timediff = state_frame.Time - state_frame.T
    grouped_state_frame = groupby(state_frame, :Time)
    day_data = combine(sdf -> sdf[argmin(sdf.Timediff), :], grouped_state_frame)
    select!(day_data, Not(:Timediff,:T))

    return day_data

end
