"""
Authors: Bennet Outland, Joshua Hoffman, McCade Hughes, Heath Buer, Joseph Romero
Project: ME 331 Final Project
Organization: South Dakota Mines
Date: 05/05/2023
License: GPLv3
"""

using CSV
using DataFrames

wt_data = CSV.File("./csv/wind_tunnel_data_frenchy.csv")

# Params
velocity = 40
α = 20

function trap_int(x, y)
    """
    Standard implementation of non-uniform trapezoid integration

    Inputs:
    * x: array of the independent variable
    * y: array of the dependent variable

    Output:
    * sum: approximate integral of y wrt x
    """

    sum = 0
    for i ∈ 2:(length(x))
        Δx = x[i] - x[i-1]
        sum = sum + ((y[i-1] + y[i])/2)*Δx
    end

    return sum
end

function find_index(wt_data, velocity, α)
    """
    Determine the indicies that correspond to the correct velocity and angle of Attack

    Inputs:
    * wt_data: wind tunnel data
    * velocity: velocity to be investigated
    * α: angle of attack to be investigated

    Output:
    * index: 
    """

    for i ∈ eachindex(wt_data["Velocity [MPH]"])
        if (wt_data["Velocity [MPH]"][i] == velocity) && (wt_data["Angle of Attack"][i] == α)
            return i
        end
    end

end

function process_data(wt_data, velocity, α; debug_msg=false)
    vel_index = find_index(wt_data, velocity, α)

    # Integrating the bottom
    bottom_indices = (10 + vel_index):(16 + vel_index)
    push!(bottom_indices, 1)

    lengths_bottom = wt_data["Distance [m]"][bottom_indices]

    bottom_extrap = wt_data["Pressure AVG"][bottom_indices[end]] * (3.5*0.0254 - lengths_bottom[end]) * 1

    bottom_force = (trap_int(lengths_bottom, wt_data["Pressure AVG"][bottom_indices]) * 1) + bottom_extrap


    # Top Integrations
    top_indices = (vel_index):(9 + vel_index)

    lengths_top = wt_data["Distance [m]"][top_indices]

    P_10 = wt_data["Pressure AVG"][top_indices[end]]
    P_9 = wt_data["Pressure AVG"][top_indices[length(lengths_top) -  1]]
    L_10 = lengths_top[end]
    L_9 = lengths_top[end - 1]
    L_end = 3.5*0.0254

    P_end = ((P_10 - P_9) / (L_10 - L_9)) * (L_end - L_10) + P_10

    top_extrap =  P_end * (L_end - L_10) + 0.5*(P_10-P_end)*(L_end - L_10)

    top_force = trap_int(lengths_top, wt_data["Pressure AVG"][top_indices]) * 1 + top_extrap

    # Calculate C_L
    ΣF = bottom_force - top_force # N
    v = velocity * 0.44704 # m/s
    A = 0.0889 # m^2
    ρ = 1.225 # kg/m^3

    C_L = (2 * ΣF) / (A * ρ * v^2)

    if debug_msg
        println("\n+===================================+")
        println("Top [N]: ", top_force)
        println("P9, P10 ", P_9, " ", P_10)
        println("Slope Test: ", ((P_10 - P_9) / (L_10 - L_9)))
        println("Bottom [N]: ", bottom_force)
        println("ΣF [N]: ", ΣF)
        println("Velocity [MPH]: ", velocity)
        println("C_L: ", C_L)
        println("α [°]: ", α)
        println("+===================================+")
    end

    return C_L, α
end

process_data(wt_data, velocity, α, debug_msg=true)