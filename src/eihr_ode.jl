"""
eihr_ode!(du, u, p, t)

Calculate the ordinary differential equations (ODEs) for the EIHR model.

Parameters:
- `du`: Array{Float64,1} - The derivative of the state variables.
- `u`: Array{Float64,1} - The current state variables.
- `p`: Tuple{Float64,Float64,Float64,Float64,Float64} - The model parameters (alpha, gamma, nu, w, epsilon).
- `t`: Float64 - The current time.

Returns:
- `du`: Array{Float64,1} - The derivative of the state variables.

"""
function eihr_ode!(du, u, p, t)
    # Prevalence ODE
    (E, I, H) = u
    (gamma, nu, epsilon, alphas, ws, param_change_times) = p

    # Time varying
    ind_t = max(searchsortedlast(param_change_times, t), 1)
    alpha = alphas[ind_t]
    w = ws[ind_t]

    # -> E
    exposed_in = alpha * I
    # E -> I
    progression = gamma * E
    # I -> H
    hospitalization = nu * w * I
    # I -> R
    non_hospitalized_recovery = nu * (1 - w) * I
    # H -> R
    hospitalized_recovery = epsilon * H

    @inbounds begin
        du[1] = exposed_in - progression # E
        du[2] = progression - (hospitalization + non_hospitalized_recovery) # I
        du[3] = hospitalization - hospitalized_recovery # H
        #du[4] = hospitalization # C_h
    end
end


function eihr_ode_inc!(du, u, p, t)
    # Incidence ODE
    (E, I, N_ih) = u
    (gamma, nu, alphas, ws, param_change_times) = p

    # Time varying
    ind_t = max(searchsortedlast(param_change_times, t), 1)
    alpha = alphas[ind_t]
    w = ws[ind_t]

    # -> E
    exposed_in = alpha * I
    # E -> I
    progression = gamma * E
    # I -> H
    hospitalization = nu * w * I
    # I -> R
    non_hospitalized_recovery = nu * (1 - w) * I
    # H -> R
    #hospitalized_recovery = epsilon * H

    @inbounds begin
        du[1] = exposed_in - progression # E
        du[2] = progression - (hospitalization + non_hospitalized_recovery) # I
        #du[3] = hospitalization - hospitalized_recovery # H
        du[3] = hospitalization # N_ih
    end
end



function eihr_ode_sim!(du, u, p, t)
    # Simulation ODE
    (E, I, H) = u
    (gamma, nu, epsilon, alphas, ws, param_change_times) = p

    # Time varying
    ind_t = max(searchsortedlast(param_change_times, t), 1)
    alpha = alphas[ind_t]
    w = ws[ind_t]

    # -> E
    exposed_in = alpha * I
    # E -> I
    progression = gamma * E
    # I -> H
    hospitalization = nu * w * I
    # I -> R
    non_hospitalized_recovery = nu * (1 - w) * I
    # H -> R
    hospitalized_recovery = epsilon * H

    @inbounds begin
        du[1] = exposed_in - progression # E
        du[2] = progression - (hospitalization + non_hospitalized_recovery) # I
        du[3] = hospitalization - hospitalized_recovery # H
        du[4] = hospitalization # C_h
    end
end

"""
eihr_ode_const_w!(du, u, p, t)

Calculate the ordinary differential equations (ODEs) for the EIHR model. With constant hospitalization probability.

Parameters:
- `du`: Array{Float64,1} - The derivative of the state variables.
- `u`: Array{Float64,1} - The current state variables.
- `p`: Tuple{Float64,Float64,Float64,Float64,Float64} - The model parameters (alpha, gamma, nu, w, epsilon).
- `t`: Float64 - The current time.

Returns:
- `du`: Array{Float64,1} - The derivative of the state variables.

"""
function eihr_ode_const_w!(du, u, p, t)
    # ODE for EIHR model with constant beta
    (E, I, H) = u
    (gamma, nu, epsilon, alphas, w, param_change_times) = p

    # Time varying
    ind_t = max(searchsortedlast(param_change_times, t), 1)
    alpha = alphas[ind_t]

    # -> E
    exposed_in = alpha * I
    # E -> I
    progression = gamma * E
    # I -> H
    hospitalization = nu * w * I
    # I -> R
    non_hospitalized_recovery = nu * (1 - w) * I
    # H -> R
    hospitalized_recovery = epsilon * H

    @inbounds begin
        du[1] = exposed_in - progression # E
        du[2] = progression - (hospitalization + non_hospitalized_recovery) # I
        du[3] = hospitalization - hospitalized_recovery # H
        du[4] = hospitalization # C_H
    end
end
