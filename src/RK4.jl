function rk3(func, t, y, dt)
    k_1 = func(t, y)
    k_2 = func(t + 0.5*dt, y + 0.5*dt*k_1)
    k_3 = func(t + 0.5*dt, y + 0.5*dt*k_2)
    k_4 = func(t + h, y + h*k_3)

    y = y + (k_1+k_4)/6 + (k_2 + k_3)/3 
    t = t+dt
end

function f()
    m = 45000 # kg
    ρ = 1.225 # kg/m^3
    A = 153.125 # m^2
    C_1 = 0.2934237  # Coefficient of Drag at 40mph
    C_2 = 0.60789517  # Coefficient of Drag at 65 mph
    V_1 = mph_2_mps(40) 
    V_2 = mph_2_mps(65) 
    δ = (C_2 - C_1)/(V_2 - V_1)   # Slope
    g = 9.81  # Gravity
    B = 4200/m; # Breaking Force

    v = u[1]
    dv = u[2]
    du[1] = dv
    du[2] =  -(ρ*A)/(2*m) * (δ*v^3 - (δ*V_1 + C_1)*v^2) - B
    return v, dv
end