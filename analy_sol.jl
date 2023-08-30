function u(x,y,z)
    a = [sin(π*i + π*j + π*k) for i ∈ x, j ∈ y, k ∈ z]
    return a
end

function u_new(x,y,z)
    a = [sin(π*i + 2π*j + 3π*k) for i ∈ x, j ∈ y, k ∈ z]
    return a
end

# x_ex = 0:0.25:1
# y_ex = 0:0.25:1
# z_ex = 0:0.25:1

# ex_analy_sol = u(x_ex,y_ex,z_ex)


function form_analy_sol(;u=u,N = 2^3)
    x_ex = 0:1/N:1
    y_ex = 0:1/N:1
    z_ex = 0:1/N:1
    ex_analy_sol = u(x_ex,y_ex,z_ex)
    return (ex_analy_sol,x_ex,y_ex,z_ex)
end


# Dirichlet conditions
function u_Front(y,z)
    a = [sin(π + π*j + π*k) for j ∈ y, k ∈ z]
    return a
end

function u_End(y,z)
    a = [sin(π*j + π*k) for j ∈ y, k ∈ z]
    return a
end

function u_Left(x,z)
    a = [sin(π*i + π*k) for i ∈ x, k ∈ z]
    return a
end

function u_Right(x,z)
    a = [sin(π*i + π*k) for i ∈ x, k ∈ z]
    return a
end

function u_Top(x,y)
    a = [sin(π*i + π*j + π) for i ∈ x, j ∈ y]
    return a
end


function u_Bottom(x,y)
    a = [sin(π*i + π*j) for i ∈ x, j ∈ y]
    return a
end


# Neumann conditions

function u_x_Front(y,z)
    a = π * [cos(π + π*j + π*k) for j ∈ y, k ∈ z]
    return a 
end


function u_x_End(y,z)
    a = -π * [cos(π*j + π*k) for j ∈ y, k ∈ z] # -1 normal direction
    return a 
end

function u_y_End(y,z)
    a = -π * [cos(π*j + π*k) for j ∈ y, k ∈ z] # -1 normal direction
    return a 
end

function u_y_Left(x,z)
    a = -π * [cos(π*i + π*k) for i ∈ x, k ∈ z] # -1 normal direction
    return a 
end

function u_x_Left(x,z)
    a = -π * [cos(π*i + π*k) for i ∈ x, k ∈ z]
end

function u_y_Right(x,z)
    a = π * [cos(π*i + π + π*k) for i ∈ x, k ∈ z] # -1 normal direction
    return a 
end

function u_x_Right(x,z)
    a = π * [cos(π*i + π + π*k) for i ∈ y, k ∈ z] # -1 normal direction
    return a 
end

function u_z_Top(x,y)
    a = π * [cos(π*i + π*j + π) for i ∈ x, j ∈ y] # -1 normal direction
    return a 
end

function u_x_Top(x,y)
    a = π * [cos(π*i + π*j + π) for i ∈ x, j ∈ y] # -1 normal direction
    return a 
end


function u_z_Bottom(x,y)
    a = -π * [cos(π*i + π*j) for i ∈ x, j ∈ y] # -1 normal direction
    return a 
end

function u_x_Bottom(x,y)
    a = -π * [cos(π*i + π*j) for i ∈ x, j ∈ y] # -1 normal direction
    return a 
end