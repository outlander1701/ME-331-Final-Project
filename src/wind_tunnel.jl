"""
Authors: Bennet Outland, Joshua Hoffman, McCade Hughes, Heath Buer, Joseph Romero
Project: ME 331 Final Project
Organization: South Dakota Mines
Date: 05/05/2023
License: GPLv3
"""

using CSV
using DataFrames
using Plots
using LaTeXStrings
include("./utilities.jl")

wt_data = CSV.File("./csv/wind_tunnel_data_frenchy.csv")
cfd_data = CSV.File("./csv/CFDForceData.csv")

# Params
velocity = 40
α = 0

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
    bottom_indices = [1]
    l_vec = (10 + vel_index):(16 + vel_index)
    for i ∈ l_vec 
        push!(bottom_indices, i)
    end

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
    bottom_force_proj = bottom_force * cos(α * pi/180)
    top_force_proj = top_force * cos(α * pi/180)

    ΣF = bottom_force_proj - top_force_proj # N
    v = velocity * 0.44704 # m/s
    A = 0.0889 * cos(α * pi/180) # m^2
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

#process_data(wt_data, velocity, α, debug_msg=true)

α_vec = [0, 5, 10, 15, 20, 0, 5, 10, 15, 20, 0, 5, 10, 15, 20, 0, 5, 10, 15, 20, 0, 5, 10, 15, 20, 0, 5, 10, 15, 20, 0, 5, 10, 15, 20]
v_vec = [40, 40, 40, 40, 40, 45, 45, 45, 45, 45, 50, 50, 50, 50, 50, 55, 55, 55, 55, 55, 60, 60, 60, 60, 60, 65, 65, 65, 65, 65, 70, 70, 70, 70, 70]
C_L_vec = []

for i ∈ eachindex(α_vec)
    push!(C_L_vec, process_data(wt_data, v_vec[i], α_vec[i], debug_msg=false)[1])
end

function plot_data_cfd_vs_α(α_vec, C_L_vec, v_vec)

    Re_dict = Dict([(40, "1.09E+05"), (45, "1.22E+05"), (50, "1.36E+05"), (55, "1.50E+05"), (60, "1.63E+05"), (65, "1.77E+05"), (70, "1.90E+05"), (90, "2.45E+05"), (200, "5.44E+05")])

    display(plot(α_vec[1:5], C_L_vec[1:5], label="Re: 1.09E+05", c=0))
    scatter!(α_vec[1:5], C_L_vec[1:5], label=nothing, c=0)
    for i ∈ 1:Int64((length(C_L_vec)/5)-1)
        velocity_label = v_vec[5*i+1]
        Re_label = Re_dict[v_vec[5*i+1]]
        if velocity_label == 40 || velocity_label == 65
            plot!(α_vec[1:5], C_L_vec[5*i+1:i*5 + 5], label="Re: $Re_label", c=i)
            scatter!(α_vec[1:5], C_L_vec[5*i+1:i*5 + 5], label=nothing, c=i)
        end
    end
   

   
    F_d = cfd_data["F_D"]
    F_l = cfd_data["F_L"]

    V = cfd_data["Velocity"]
    α_cfd = cfd_data["AoA"]

    C_d = Coefficient(F_d, V)
    C_l = Coefficient(F_l, V)

    Cd_sorted = reshape(C_d, 4, 8)
    Cl_sorted = reshape(C_l, 4, 8)
    α_sorted = reshape(α_cfd, 4, 8)

    N = length(α_sorted[:,1])

    for i ∈ 1:N
        velocity_label = V[i]
        Re_label = Re_dict[V[i]]
        if velocity_label == 40 || velocity_label == 65
            plot!(α_sorted[i,1:8], Cl_sorted[i,1:8], label="Re: $Re_label (CFD)", c=i)
            scatter!(α_sorted[i,1:8], Cl_sorted[i,1:8], label=nothing, c=i)
        end
    end
  

    xlabel!(L"Angle \, \, of \, \, Attack \, \,  (α) \, \, [deg]")
    ylabel!(L"Coefficent \, \, of \, \, Lift \, \, (C_{L})")
end

function plot_data_cl_α(cfd_data)

    Re_dict = Dict([(40, "1.09E+05"), (45, "1.22E+05"), (50, "1.36E+05"), (55, "1.50E+05"), (60, "1.63E+05"), (65, "1.77E+05"), (70, "1.90E+05"), (90, "2.45E+05"), (200, "5.44E+05")])

    F_d = cfd_data["F_D"]
    F_l = cfd_data["F_L"]

    V = cfd_data["Velocity"]
    α_cfd = cfd_data["AoA"]

    C_d = Coefficient(F_d, V)
    C_l = Coefficient(F_l, V)

    Cd_sorted = reshape(C_d, 4, 8)
    Cl_sorted = reshape(C_l, 4, 8)
    α_sorted = reshape(α_cfd, 4, 8)

    N = length(α_sorted[:,1])

    plot()
    for i ∈ 1:N
        velocity_label = V[i]
        Re_label = Re_dict[V[i]]

        plot!(α_sorted[i,1:8], Cl_sorted[i,1:8], label="Re: $Re_label")
        scatter!(α_sorted[i,1:8], Cl_sorted[i,1:8], label=nothing, c=i)
        
    end


    xlabel!(L"Angle \, \, of \, \, Attack \, \,  (α) \, \, [deg]")
    ylabel!(L"Coefficent \, \, of \, \, Lift \, \, (C_{L})")
end

function plot_data_cd_α(cfd_data)

    Re_dict = Dict([(40, "1.09E+05"), (45, "1.22E+05"), (50, "1.36E+05"), (55, "1.50E+05"), (60, "1.63E+05"), (65, "1.77E+05"), (70, "1.90E+05"), (90, "2.45E+05"), (200, "5.44E+05")])

    F_d = cfd_data["F_D"]
    F_l = cfd_data["F_L"]

    V = cfd_data["Velocity"]
    α_cfd = cfd_data["AoA"]

    C_d = Coefficient(F_d, V)
    C_l = Coefficient(F_l, V)

    Cd_sorted = reshape(C_d, 4, 8)
    Cl_sorted = reshape(C_l, 4, 8)
    α_sorted = reshape(α_cfd, 4, 8)

    N = length(α_sorted[:,1])

    plot()
    for i ∈ 1:N
        velocity_label = V[i]
        Re_label = Re_dict[V[i]]

        plot!(α_sorted[i,1:8], Cd_sorted[i,1:8], label="Re: $Re_label", c=i)
        scatter!(α_sorted[i,1:8], Cd_sorted[i,1:8], label=nothing, c=i)
        
    end


    xlabel!(L"Angle \, \, of \, \, Attack \, \,  (α) \, \, [deg]")
    ylabel!(L"Coefficent \, \, of \, \, Drag \, \, (C_{D})")
end

function plot_data_cl_cd_α(cfd_data)

    Re_dict = Dict([(40, "1.09E+05"), (45, "1.22E+05"), (50, "1.36E+05"), (55, "1.50E+05"), (60, "1.63E+05"), (65, "1.77E+05"), (70, "1.90E+05"), (90, "2.45E+05"), (200, "5.44E+05")])

    F_d = cfd_data["F_D"]
    F_l = cfd_data["F_L"]

    V = cfd_data["Velocity"]
    α_cfd = cfd_data["AoA"]

    C_d = Coefficient(F_d, V)
    C_l = Coefficient(F_l, V)

    Cd_sorted = reshape(C_d, 4, 8)
    Cl_sorted = reshape(C_l, 4, 8)
    α_sorted = reshape(α_cfd, 4, 8)

    ratio_sorted = Cl_sorted ./ Cd_sorted

    N = length(α_sorted[:,1])

    plot()
    for i ∈ 1:N
        velocity_label = V[i]
        Re_label = Re_dict[V[i]]

        plot!(α_sorted[i,1:8], ratio_sorted[i,1:8], label="Re: $Re_label", c=i)
        scatter!(α_sorted[i,1:8], ratio_sorted[i,1:8], label=nothing, c=i)
        
    end


    xlabel!(L"Angle \, \, of \, \, Attack \, \,  (α) \, \, [deg]")
    ylabel!(L"(C_{L} / C_{D})")
end

function plot_data_cl_cd(cfd_data)

    Re_dict = Dict([(40, "1.09E+05"), (45, "1.22E+05"), (50, "1.36E+05"), (55, "1.50E+05"), (60, "1.63E+05"), (65, "1.77E+05"), (70, "1.90E+05"), (90, "2.45E+05"), (200, "5.44E+05")])

    F_d = cfd_data["F_D"]
    F_l = cfd_data["F_L"]

    V = cfd_data["Velocity"]
    α_cfd = cfd_data["AoA"]

    C_d = Coefficient(F_d, V)
    C_l = Coefficient(F_l, V)

    Cd_sorted = reshape(C_d, 4, 8)
    Cl_sorted = reshape(C_l, 4, 8)
    α_sorted = reshape(α_cfd, 4, 8)

    ratio_sorted = Cl_sorted ./ Cd_sorted

    N = length(α_sorted[:,1])

    plot()
    for i ∈ 1:N
        velocity_label = V[i]
        Re_label = Re_dict[V[i]]

        plot!(Cd_sorted[i,1:8], Cl_sorted[i,1:8], label="Re: $Re_label", c=i)
        scatter!(Cd_sorted[i,1:8], Cl_sorted[i,1:8], label=nothing, c=i)
        
    end


    xlabel!(L"C_{D}")
    ylabel!(L"C_{L}")
end

#plot_data_cfd_vs_α(α_vec, C_L_vec, v_vec)
#plot_data_cl_α(cfd_data)
#plot_data_cd_α(cfd_data)
#plot_data_cl_cd_α(cfd_data)
plot_data_cl_cd(cfd_data)