(* ::Package:: *)

(*             GENERAL INPUT              *)
(* number of bath sites *)
Nbath = 3; 
(* number of orbitals *)
Norb = 1; 
(* number of impurity sites *)
Nimp = 1; 
(* total number of sites: bath+impurity *)
L = Nimp + Nbath; 
(* number of spin states *)
f = 2; 
(* call the function EdModeInfo[EdMode] to get details *)
EdMode = "Raman"; 
(* lattice crystal structure: "Bethe", "Hypercubic", etc. *)
LatticeType = "Hypercubic"; 
(* lattice dimensionality *)
LatticeDim = 1; 
(* lattice points *)
LatticePoints = 1000;
(* compute many-body functions independently for two different sublattices A and B? *)
SublatticesQ = False;
(* half of the on-site energy split between the two sublattices *)
V = 0.0;
(* set True to enforce orbital symmetry and avoid repeating calculations *)
OrbitalSymmetry = True; 
(* list of half-bandwidths for all the orbitals *)
W = ConstantArray[1.0, Norb];

(* add orbital hybridization options if needed ... *)


(*      INPUT PHYSICAL PARAMETERS        *)
(* interaction energy (list of U values for the orbitals) *)
U = ConstantArray[0.0, Norb]; 
(* Hund's J. It's used only when HundMode = True to enforce rotation invariance of the Kanamori model. *)
JH = 0.0; 
(* density-density opposite spin coupling. It is set automatically if HundMode = True. *)
Ust = 0.00; 
(* density-density same spin coupling. It is set automatically if HundMode = True. *)
Usec = 0.00; 
(* pair-hopping coupling. It is set automatically if HundMode = True. *)
Jph = 0.0 * IdentityMatrix[Norb]; 
(* spin-exchange coupling. It is set automatically if HundMode = True. *)
Jse = 0.0; 
(* if this is True, interorbital couplings are authomatically set to Ust=U-2JH; Usec=U-3JH, Jph=JH, Jse=-JH enforcing the rotational invariance *)
HundMode = False; 
(* chemical potential *)
\[Mu] = 0.0; 
(* Crystal field splitting *)
\[Delta] = {0};
(* Explicit magnetic field: (can be different for different orbitals) *)
h = ConstantArray[0.0, {Norb, f}];
(* Raman hopping matrix (can be different for different orbitals) *)
M = 2.5 * ConstantArray[PauliMatrix[1], Norb];
(* Gauge field associated with Raman tunneling *)
\[Gamma] = 0.00 * Pi;
(* unit vector giving the spatial direction where tunneling is modified by the gauge field *)
u = {1.0}; 
(* temperature *)
T = 0; 


(* INFO ON MATSUBARA AND REAL FREQUENCIES *)
(* fictitious temperature to define Matsubara frequencies *)
TMats = 0.001; 
(* Total number of Matsubara frequencies *)
NMatsubara = 5000; 
(* set min and max value for the set of real frequencies *)
\[Omega]min = -5.; \[Omega]max = 5.; 
(* number of real frequencies *)
NReal = 10000;
(* real frequency step *)
d\[Omega] = (\[Omega]max - \[Omega]min)/NReal; 
(* small shift of the pole in the imaginary axis: this avoids singularities, but introduces an artificial broadening of the spectrum *)
\[Eta] = 0.02;


(* INPUT-OUTPUT MANAGEMENT *)
(* directory where the code is *)
CodeDirectory = NotebookDirectory[];
(* directory where the output will be stored *)
OutputDirectory = CodeDirectory<>"prova\\";(*"Jph="<>ToString[Jph[[1,2]]]<>"\\";*)
(* load Hamiltonian from a file? *)
LoadHamiltonianQ = False;


(* NUMERICAL DETAILS OF THE ALGORITHM *)
(* path to input file of bath parameters or "Default" *)
InitializeBathMode = "Default";
(* if this is True, an extra on-site energy term is added to ensure PH symmetry of the Kanamori Hamiltonian when \[Mu]=0 *)
HFMode = True; 
(* use with caution: if this is True, the full spectrum will be automatically determined! *)
FullDiagonalizationMode = True;
(* below this threshold, two energy levels are assumed to be degenerate *)
DegeneracyThreshold = 10^(-8);
(* dimension threshold above which Lanczos is used to perform exact diagonalization *)
MinLanczosDim = 2000;
(* maximum number of Lanczos iterations for spectrum determination *)
MaxLanczosIter = 2000;
(* if you use Lanczos for spectrum determination, compute AT LEAST this number of eigenstates *)
MinNumberOfEigs = 10;
(* minimum number of DMFT loops *)
DMFTMinIterations = 2; 
(* maximum number of DMFT loops *)
DMFTMaxIterations = 30; 
(* threshold for DMFT loop convergence *)
DMFTerror = 1.0 * 10^(-5); 
(* InverseG0 = Mixing * InverseG0old + (1 - Mixing) * InverseG0 *)
Mixing = 0.25; 
(* type of minimization: "Global" or "Local" *)
MinimizationType = "Local";
(* Method for the minimization procedure *)
MinimizationMethod = "ConjugateGradient";
(* Max number of Conjugate Gradient iterations *)
CGMaxIterations = 1000;
(* Number of Matsubara frequencies used to perform the fit *)
CGNMatsubara = 1500;
(* Accuracy goal for minimization procedure *)
CGAccuracy = 7;
(* list of weights attributed to each Matsubara frequency w(i\[Omega]): IT SHOULD BE >= CGNMatsubara *)
CGWeight = ConstantArray[1., CGNMatsubara];

(* turn off annoying messages from Eigensystem and stop evaluation after the first error occurs *)
AbortAtFirstErrorQ = False;


Print["Input file imported successfully. "]
