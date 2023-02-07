using Symbolics

@variables x y z

π_s = Symbolics.pi

u1 = cos(2*π_s*x) * cos(π_s*y) * cos(π_s*z)
u2 = cos(π_s*x) * cos(3*π_s*y) * cos(π_s*z)
u3 = cos(π_s*x) * cos(π_s*y) * cos(π_s*4z)

Dx = Differential(x)
Dy = Differential(y)
Dz = Differential(z)
Dxx = Differential(x)^2
Dyy = Differential(y)^2
Dzz = Differential(z)^2
Dxy = Dx * Dy
Dxz = Dx * Dz
Dyz = Dy * Dz
Dzx = Dxz
Dzy = Dyz
Dyx = Dxy

expand_derivatives(Dx(u1)) / π_s
expand_derivatives(Dy(u1)) / π_s

expand_derivatives(Dxy(u1)) / π_s^2
expand_derivatives(Dyx(u1)) / π_s^2
expand_derivatives(Dxy(u2)) / π_s^2

expand_derivatives(Dxx(u1)) / π_s^2