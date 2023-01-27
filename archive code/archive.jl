# H and HI Operators for End and Front
# Size: N^3 by N^3
HI_Front = kron(I_Nz,I_Ny,HIx)
HI_End = kron(I_Nz,I_Ny,HIx)
H_Front = kron(I_Nz, I_Ny, H1x)
H_End = kron(I_Nz, I_Ny, H1x)

# H and HI operators for Left and Right
# Size: N^3 by N^3
HI_Left = kron(I_Nz,HIy,I_Nx)
HI_Right = kron(I_Nz,HIy,I_Nx)
H_Left = kron(I_Nz, H1y, I_Nx)
H_Right = kron(I_Nz, H1y, I_Nx)

# H and HI operators for Bottom and Top
# Size : N^3 by N^3
HI_Bottom = kron(HIz,I_Nx,I_Ny)
HI_Top = kron(HIz,I_Nx,I_Ny)
H_Bottom = kron(H1z, I_Ny, I_Nx)
H_Top = kron(H1z, I_Ny, I_Nx)
