using DifferentialEquations
using Plots

# Guide: https://docs.sciml.ai/DiffEqDocs/stable/examples/classical_physics/#Second-Order-Non-linear-ODE

# This is still very much WIP do not run

# Constants
m = 42
ρ = 1000
A = 42
C_D = 42
C_L = 42
μ = 42
g = 9.81
μ = 0.3


# Initial Conditions
u₀ = [0, π / 2]
t_span = (0.0, 120)

# Define the problem
function airplane_landing(du, u, p, t)
    """
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)
    """

    x_dot = u[1]
    
    
    
end

#Pass to solvers
prob = ODEProblem(airplane_landing, u₀, t_span)
sol = solve(prob, Tsit5())

#Plot
plot(sol, linewidth = 2, title = "Landing Problem", xaxis = "Time",
     yaxis = "Distance", label = ["x" "dx"])
