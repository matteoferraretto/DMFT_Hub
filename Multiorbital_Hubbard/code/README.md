The code is organized as follows. 
There are 10 packages:

  - MyLinearAlgebra`
  - Parameters`
  - Sectors`
  - Hamiltonian`
  - Observables`
  - ImpurityGreenFunction`
  - Lattices`
  - SelfConsistency`
  - Facilities`
  - DMFT`

These packages are essentially libraries where the functions are stored, with the exception of DMFT`, which is just used as a setup code.

There are other 5 files:

  - InputFile_Template.wl
  - Preparation.wl
  - DMFT_Loop.wl
  - Post_Processing.wl
  - Main.wl

The file InputFile_Template.wl contains a list of tunable input variables, so the first step to run the code is to modify this file according to the users' needs.
The file Preparation.wl is used to initialize sectors, hamiltonians, lattices and bath parameters and to perform relevant checks on the input consistency. 
The file DMFT_Loop.wl performs the DMFT iterations
The file Post_Processing.wl is called after convergence to plot and print information obtained with a standardized data analysis. 
The file Main.wl calls the other three in order and displays the output of the caltulation. To launch the code you just need to open Main.wl and run all the cells with the command "Evaluate Notebook".
