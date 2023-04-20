using CSV
include("utilities.jl")

# Import 
cfd_data = CSV.File("./csv/cfd_data.csv")

vel_array = []
C_L_vec = []

for i ∈ eachindex(cfd_data["Angle"])
    if cfd_data["Angle"][i] == 0
        push!(vel_array, mph_2_mps(cfd_data["Velocity"][i]))
        push!(C_L_vec, cfd_data["Cl"][i])
    end
end

function find_val(vel_array)
    
    ϵ = 0.1
    for j ∈ minimum(vel_array):ϵ:maximum(vel_array)
        C_L_est = lin_interp(j, vel_array, C_L_vec)
        vel_est = predicted_vel(C_L_est)
        if abs(vel_est - j) < 1
            return [vel_est, C_L_est]
        end
    end

end

println(find_val(vel_array))