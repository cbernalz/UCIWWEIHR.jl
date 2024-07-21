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
    # ODE for EIHR model with constant beta
    (E, I, H) = u
    (alpha, gamma, nu, w, epsilon) = p
    # w - ratio of infected individuals going to hospital

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
    end
end