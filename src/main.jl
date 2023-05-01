using CSV
include("utilities.jl")
include("CFD.jl")

# Constants
m = 45000 # kg
g = 9.81 # m/s^2
A = 153.125 # m^2
ρ = 1.225 # kg/m^3

# Import 
cfd_data = CSV.File("./csv/CFDForceData.csv")

vel_array = V[1:4] .* 0.44704

println("+=========================================+")
println("+============== Input Check ==============+")
println("+=========================================+")
println(" ")

println("Velocities: ", vel_array)
println("C_L: ", Cl_sorted[:,3])
println("C_D: ", Cd_sorted[:,3])

println(" ")

C_L_vec = []
C_D_vec = []

function find_vel(vel_array, C_L, C_D)
    C_L_vec = C_L[:,3]
    C_D_vec = C_D[:,3]

    ϵ = 0.001
    for j ∈ minimum(vel_array):ϵ:maximum(vel_array)
        C_L_est = lin_interp(j, vel_array, C_L_vec)
        vel_est = predicted_vel(C_L_est)

        if abs(vel_est - j) < 1
            C_D_est = lin_interp(vel_est, vel_array, C_D_vec)
            return vel_est, C_L_est, C_D_est
        end
    end

end


v, C_L, C_D = find_vel(vel_array, Cl_sorted, Cd_sorted)

F_D = 0.5*C_D*ρ*(v^2)*A

t = (m*v)/F_D

Δx = v*t + (F_D/m)*(t^2); #x = m / (C_D * ρ * A)

a = (v^2)/(2*Δx)

F_t = m*a + F_D #F_t = ((m * v^2)/(2 * x)) + ((C_D * ρ * v^2 * A) / 2)

println("+======================================+")
println("+=========== Output Values ============+")
println("+======================================+")
println(" ")

println("C_L = ", C_L)
println("C_D = ", C_D)
println("V [m/s] = ", v)
println("V [mph] = ", v/0.44704)
println("Landing Time [s] = ", t)
println("Runway length [m]: ", Δx)
println("Thrust [kN]: ", F_t/1000)

