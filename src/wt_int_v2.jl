using CSV
using DataFrames

wt_data = CSV.File("./csv/wind_tunnel_test.csv")

# Params
velocity = 50
α = 0

function trap_int(x, y)
    sum = 0
    for i ∈ 2:(length(x))
        Δx = x[i] - x[i-1]
        sum = sum + ((y[i-1] + y[i])/2)*Δx
    end

    return sum
end

# Find the data
function find_index(wt_data, velocity, α)
    
    for i ∈ eachindex(wt_data["Velocity [MPH]"])
        if (wt_data["Velocity [MPH]"][i] == velocity) && (wt_data["Angle of Attack"][i] == α)
            return i
        end
    end

end

vel_index = find_index(wt_data, velocity, α)


# Integrating the bottom
bottom_indices = (10 + vel_index):(16 + vel_index)

lengths_bottom = wt_data["Distance [in]"][bottom_indices]

bottom_extrap = wt_data["Pressure AVG [psi]"][bottom_indices[end]] * (3.5*0.0254 - lengths_bottom[end]) * 1

bottom_force = (trap_int(lengths_bottom, wt_data["Pressure AVG [psi]"][bottom_indices]) * 1) + bottom_extrap


# Top Integrations
top_indices = (1 + vel_index):(9 + vel_index)

lengths_top = wt_data["Distance [in]"][top_indices]

P_10 = wt_data["Pressure AVG [psi]"][top_indices[end]]
P_9 = wt_data["Pressure AVG [psi]"][top_indices[length(lengths_top) -  1]]
L_10 = lengths_top[end]
L_9 = lengths_top[end - 1]
L_end = 3.5*0.0254

println("Pressure ", P_9, " ", P_10)
println("Slope Test: ", ((P_10 - P_9) / (L_10 - L_9)))

P_end = ((P_10 - P_9) / (L_10 - L_9)) * (L_end - L_10) + P_10

top_extrap =  P_end * (L_end - L_10) + 0.5*(P_10-P_end)*(L_end - L_10)

top_force = trap_int(lengths_top, wt_data["Pressure AVG [psi]"][top_indices]) * 1 + top_extrap

# Calculate C_L
ΣF_freedom = bottom_force - top_force # lbf
#ΣF_frenchy = ΣF_freedom * 4.448222 # N
v = velocity * 0.44704 # m/s
A = 0.0889 # m^2
ρ = 1.225 # kg/m^3

C_L = (2 * ΣF_freedom) / (A * ρ * v^2)

println("\n+===================================+")
println("Top [lbf]: ", top_force)
println("Bottom [lbf]: ", bottom_force)
println("ΣF [N]: ", ΣF_freedom)
println("Velocity [MPH]: ", velocity)
println("α [°]: ", α)
println("C_L: ", C_L)
