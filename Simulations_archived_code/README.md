# Simulations Folder

This folder attemps to build code for BP5 simulation https://strike.scec.org/cvws/seas/download/SEAS_BP5-QD-FD.pdf
or [BP5 problem description](../Notebooks/SEAS_BP5-QD-FD.pdf)

- ``Assembling_3D_matrices.jl`` contains code to assemble linear system for the motion governed by momentum balance. The boundary conditions and notations are included in ``"../3D_face.jl"`` file.
- ``helper.jl`` contains helper functions for calculations used in SEAS BP5 benchmark problem
- ``odefun.jl`` contains formulation for the odefunc to pass to the ``ODEProblem()`` in ``DifferentialEquations.jl``. I am still designing params and odefunc to work with the ``ODEProblem()`` interface.
- ``domain.jl`` contains information to form the computational domain. Right now it's just a 256 by 256 by 256 grid in 3D
- ``coefficients.jl`` contains coefficients used in BP5 problem. The values are included in the problem description
- ``BP5-QD.jl`` contains the main function for the simulation. It is the file that should be run when doing the simulation