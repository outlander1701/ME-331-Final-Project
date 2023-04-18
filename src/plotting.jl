using CSV

# Import 
cfd_data = CSV.read("./csv/cfd.csv")
wt_data = CSV.read("./csv/wind_tunnel.csv")

data = {"CFD" => cfd_data, "Wind Tunnel" => wt_data}

plot_type = "C_L" # options: "C_L", "C_D", "C_L/C_D"

if plot_type == "C_L"
    plot(data["CFD"]["alpha"], data["CFD"]["C_L"])
elseif plot_type == "C_D"
    plot(data["CFD"]["alpha"], data["CFD"]["C_D"])
elseif plot_type == "C_L/C_D"
    plot(data["CFD"]["alpha"], data["CFD"]["C_L"] / data["CFD"]["C_D"])
else
    throw(MethodError("That is not a valid option. Options include: \"C_L\", \"C_D\", \"C_L/C_D\""))
