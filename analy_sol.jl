function u(x,y,z)
    a = [sin(π*i + π*j + π*k) for i ∈ x, j ∈ y, k ∈ z]
    return a
end

function u2(x,y,z)
    a = [sin(π*i + 2π*j + 3π*k) for i ∈ x, j ∈ y, k ∈ z]
    return a
end

x_ex = 0:0.25:1
y_ex = 0:0.25:1
z_ex = 0:0.25:1

ex_analy_sol = u(x_ex,y_ex,z_ex)


function form_analy_sol(;u=u,N = 2^3)
    x_ex = 0:1/N:1
    y_ex = 0:1/N:1
    z_ex = 0:1/N:1
    ex_analy_sol = u(x_ex,y_ex,z_ex)
    return (ex_analy_sol,x_ex,y_ex,z_ex)
end