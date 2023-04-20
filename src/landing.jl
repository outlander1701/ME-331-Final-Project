using DifferentialEquations
using Plots
include("utilities.jl")

# Guide: https://docs.sciml.ai/DiffEqDocs/stable/examples/classical_physics/#Second-Order-Non-linear-ODE

# This is still very much WIP do not run

# Constants

# Initial Conditions
v_0 = 27.662574173754745
C_D_0 = 0.5686416841355989

ρ = 1.225
A = 153.125 # m^2
m = 45000

F_D = 0.5(C_D * ρ * A * v_0^2)
F_B = 42000
a_0 = (-F_D - F_B)/m

u₀ = [v_0, a_0] # liftoff velocity, solved initial acceleration
t_span = (0., 40)

# Define the problem
function airplane_landing(du, u, p, t)
    """
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)
    """
    m = 45000 # kg
    ρ = 1.225 # kg/m^3
    A = 153.125 # m^2
    C_1 = 0.2934237  # Coefficient of Drag at 40mph
    C_2 = 0.60789517  # Coefficient of Drag at 65 mph
    V_1 = mph_2_mps(40) 
    V_2 = mph_2_mps(65) 
    δ = (C_2 - C_1)/(V_2 - V_1)   # Slope
    g = 9.81  # Gravity
    B = 4200/m; # Breaking Force
    """
    a = u[1]
    v = u[2]
    du[1] = v
    du[2] =  -(ρ*A)/(2*m) * (δ*v^3 - (δ*V_1 + C_1)*v^2) - B
    """
    v = u[1]
    dv = u[2]
    du[1] = dv
    du[2] =  -(ρ*A)/(2*m) * (δ*v^3 - (δ*V_1 + C_1)*v^2) - B
end


#Plot
#plot(sol, linewidth = 2, title = "Landing Problem", xaxis = "Time", yaxis = "Velocity", label = ["x" "dx"])

#Pass to solvers
prob = ODEProblem(airplane_landing, u₀, t_span)
sol = solve(prob, Tsit5())

println(sol.u[:,1])



#println(sol[1, 1])
dt = (t_span[2] - t_span[1]) / length(sol)
tspan = 0:dt:t_span[2]
search_index = 0;
for i ∈ eachindex(sol)
    if (sol[2, i] <= 0)
        global search_index = i-1;
        break
    end
end

#vel = sol(my_time)[1]

vel_vec = sol.(0:dt:t_span[2])

plot(sol.t, sol.u[:,1])
#function trapezoid()
