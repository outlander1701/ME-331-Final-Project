using CSV, Plots

include("./utilities.jl")

Data = CSV.File("./csv/CFDForceData.csv")

F_d = Data["F_D"]
F_l = Data["F_L"]

V = Data["Velocity"]
α = Data["AoA"]

C_d = Coefficient(F_d, V)
C_l = Coefficient(F_l, V)

Cd_sorted = reshape(C_d, 4, 8)
Cl_sorted = reshape(C_l, 4, 8)
α_sorted = reshape(α, 4, 8)

#plot_c_vs_α(α_sorted, Cl_sorted)

#plot_c_vs_α(α_sorted, Cd_sorted)

print(C_l)
#print(V)