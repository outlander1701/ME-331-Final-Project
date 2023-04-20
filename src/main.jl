using CSV
include("utilities.jl")

# Constants
m = 45000 # kg
g = 9.81 # m/s^2
A = 153.125 # m^2
ρ = 1.225 # kg/m^3

# Import 
cfd_data = CSV.File("./csv/cfd_data.csv")

vel_array = []
C_L_vec = []
C_D_vec = []

for i ∈ eachindex(cfd_data["Angle"])
    if cfd_data["Angle"][i] == 0
        push!(vel_array, mph_2_mps(cfd_data["Velocity"][i]))
        push!(C_L_vec, cfd_data["Cl"][i])
        push!(C_D_vec, cfd_data["Cd"][i])
    end
end

function find_vel(vel_array)
    
    ϵ = 0.1
    for j ∈ minimum(vel_array):ϵ:maximum(vel_array)
        C_L_est = lin_interp(j, vel_array, C_L_vec)
        vel_est = predicted_vel(C_L_est)
        if abs(vel_est - j) < 1
            C_D_est = lin_interp(vel_est, vel_array, C_D_vec)
            return vel_est, C_L_est, C_D_est
        end
    end

end

v, C_L, C_D = find_vel(vel_array)

x = m / (C_D * ρ * A)

F_t = ((m * v^2)/(2 * x)) + ((C_D * ρ * v^2 * A) / 2)


println(x, " ", F_t)