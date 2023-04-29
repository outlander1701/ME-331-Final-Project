using CSV
using DataFrames

wt_data = CSV.File("./csv/wind_tunnel_data.csv")

length_dict = Dict([["1", 0], ["2", 0.175], ["3", 0.35], ["4", 0.7], ["5", 1.05], ["6", 1.4], ["7", 1.75], ["8", 2.1], ["9", 2.45], ["10", 2.8], ["11", 0.175], ["12", 0.35], ["13", 0.7], ["14", 1.05], ["15", 1.4], ["16", 1.75], ["17", 2.1], ["18", 2.45]])

bottom = []
for i ∈ 11:17
    push!(bottom, "$i")
end

top = []
for i ∈ 2:10
    push!(top, "$i")
end

# Params
velocity = 60
α = 0
ρ = 0.0765 #(lbm/ft^3) https://www.engineersedge.com/calculators/air-density.htm#:~:text=At%20sea%20level%20and%20at,ISA%20(International%20Standard%20Atmosphere
A = 3.5/144 # ft^2

function trap_int(x, y)
    Δx = (x[end]-x[1])/length(x)

    sum = Δx*0.5*y[1] + Δx*0.5*y[end]
    #println("Init: ", sum)

    for i ∈ 2:(length(x)-1)
        sum = sum + Δx*y[i]
    end

    #println("Fin: ", sum)

    return sum

end

# Find the data

function find_index(wt_data, velocity, α)

    for i ∈ eachindex(wt_data["Velocity [MPH]"])
        if wt_data["Velocity [MPH]"][i] == velocity
            if wt_data["Angle of Attack"][i] == α
                return i
            end
        end
    end

end

function find_Cl(wt_data, velocity, α, ρ, A)
    vel_index = find_index(wt_data, velocity, α)

    # Integrating the bottom

    bottom_indices = (10 + vel_index):(17 + vel_index)

    lengths = []
    for key in bottom
    push!(lengths, length_dict["$key"])
    end

    bottom_extrap = wt_data["Pressure AVG [psi]"][bottom_indices[end]] * (3.5 - length_dict["18"]) * 1

    bottom_force = (trap_int(lengths, wt_data["Pressure AVG [psi]"][bottom_indices]) * 1) + bottom_extrap

    # Integrating the top

    top_indices = (1 + vel_index):(9 + vel_index)

    lengths = []
    for key in top
    push!(lengths, length_dict["$key"])
    end

    L_10 = length_dict["10"]
    m_top = ((wt_data["Pressure AVG [psi]"][top_indices[9]] - wt_data["Pressure AVG [psi]"][top_indices[8]]) / (L_10 - length_dict["9"]))

    top_extrap = wt_data["Pressure AVG [psi]"][top_indices[9]] * (3.5 - L_10) + 0.5 * m_top *(3.5 - L_10) * 1

    top_force = trap_int(lengths, wt_data["Pressure AVG [psi]"][top_indices]) * 1 + top_extrap

    # Calculate C_L
    ΣF = top_force - bottom_force
    v = velocity * (5280)/(60*60) # ft/s

    C_L = (2 * ΣF) / (A * ρ * v^2)

    println("\n+===================================+")
    println("ΣF [lbf]: ", ΣF)
    println("Velocity [MPH]: ", velocity)
    println("α [°]: ", α)
    println("C_L: ", C_L)
end



