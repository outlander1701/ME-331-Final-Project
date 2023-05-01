function lin_interp(x, x_vec::Vector, y_vec::Vector)
    """
    Interpolates between two points to find the y value given an x

    Inputs:

    * x: desired x value to interpolate

    * x_vec: vector of x to search through

    * y_vec: vector of y to search through

    Output:

    * y: interpolated value
    """

    n::Int8 = -42
    for i ∈ 1:length(x_vec)
        if x_vec[i] > x
            n = i
            break
        end
    end

    return (((y_vec[n] - y_vec[n-1]) / (x_vec[n] - x_vec[n-1])) * (x - x_vec[n-1])) + y_vec[n-1]

end

function predicted_vel(C_l)
    m = 45000 # kg
    g = 9.81
    A = 153.125 # m^2
    ρ = 1.225 # kg/m^3

    v2 = (2*m*g) / (C_l * A * ρ)
    
    if v2 > 0
        return sqrt(v2)
    else
        return 0
    end

end

function mps_2_mph(speed)
    return 2.236936 * speed
end

function mph_2_mps(speed)
    return (1/2.236936) * speed
end

function Coefficient(F, V)
    
    N = length(F);
    C = Array{Float16, 1}(undef, N);
    V_metric = V .* 0.44704

    for i ∈ 1:N
        C[i] = F[i]/(0.5*1.225*(V_metric[i]^2)*(0.0889))
    end
    return C
end

function plot_c_vs_α(α, Coefficient)

    N = length(α[:,1])

    plot()

    for i ∈ 1:N
        display(plot!(α[i,1:8], Coefficient[i,1:8]))
    end

end