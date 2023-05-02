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

function find_vel_landing(vel_array, C_L, C_D)
    C_L_vec = C_L[:,7]
    C_D_vec = C_D[:,7]

    ϵ = 0.001
    for j ∈ minimum(vel_array):ϵ:maximum(vel_array)
        C_L_landing = lin_interp(j, vel_array, C_L_vec)
        vel_landing = predicted_vel(C_L_landing)
        
        if abs(vel_landing - j) < 1
            C_D_landing = lin_interp(vel_landing, vel_array, C_D_vec)
            return vel_landing, C_L_landing, C_D_landing
        end
    end

end
function find_vel_takeoff(vel_array, C_L, C_D)
    C_L_vec = C_L[:,6]
    C_D_vec = C_D[:,6]
    C_L_takeoff = 0

    ϵ = 0.001

    for j ∈ minimum(vel_array):ϵ:maximum(vel_array)
        C_L_last = C_L_takeoff

        C_L_takeoff = lin_interp(j, vel_array, C_L_vec)
        vel_takeoff = predicted_vel(C_L_takeoff)
        
        if (abs(vel_takeoff - j) < 1 && C_L_last > C_L_takeoff)
            C_D_takeoff = lin_interp(vel_takeoff, vel_array, C_D_vec)
            return vel_takeoff, C_L_takeoff, C_D_takeoff
        end
    end

end

# Landing Drag force, time, and distance
v_land, C_L_land, C_D_land = find_vel(vel_array, Cl_sorted, Cd_sorted)

F_D = 0.5*C_D_land*ρ*(v_land^2)*A

t_land = (m*v_land)/F_D

Δx_land = v_land*t + (F_D/m)*(t^2); #x = m / (C_D * ρ * A)

# Takeoff acceleration, distance, thrust
v_takeoff, C_L_takeoff, C_D_takeoff = find_vel_takeoff(vel_array, Cl_sorted, Cd_sorted)

t_takeoff = 30; # [s]

F_D = 0.5*C_D_takeoff*ρ*(v_takeoff^2)*A

a_takeoff = v_takeoff/t #a = (v_takeoff^2)/(2*Δx_land)

Δx_takeoff = 0.5*a_takeoff*(t_takeoff^2)

F_t = m*a + F_D #F_t = ((m * v^2)/(2 * x)) + ((C_D * ρ * v^2 * A) / 2)

println("+======================================+")
println("+============== Landing ===============+")
println("+======================================+")
println(" ")

println("C_L = ", C_L_land)
println("C_D = ", C_D_land)
println("V [m/s] = ", v_land)
println("V [mph] = ", v_land/0.44704)
println("Landing Time [s] = ", t_land)
println("Length Required [m]: ", Δx_land)

println(" ")
println("+======================================+")
println("+============== Takeoff ===============+")
println("+======================================+")
println(" ")

println("C_L = ", C_L_takeoff)
println("C_D = ", C_D_takeoff)
println("V [m/s] = ", v_takeoff)
println("V [mph] = ", v_takeoff/0.44704)
println("Takeoff Time [s] = ", t_takeoff)
println("Length Required [m]: ", Δx_takeoff)
println("Thrust Required [kN]: ", F_t/1000)

