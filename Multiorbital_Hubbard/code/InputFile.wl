(* ::Package:: *)

(*             GENERAL INPUT              *)
(* number of bath sites *)
Nbath = 4; 
(* number of orbitals *)
Norb = 2; 
(* number of impurity sites *)
Nimp = 1; 
(* total number of sites: bath+impurity *)
L = Nimp + Nbath; 
(* number of spin states *)
f = 2; 
(* call the function EdModeInfo[EdMode] to get details *)
EdMode = "Magnetic"; 
(* lattice crystal structure: "Bethe", "Hypercubic", etc. *)
LatticeType = "Hypercubic"; 
(* lattice dimensionality *)
LatticeDim = 2; 
(* lattice points *)
LatticePoints = 1600;
(* compute many-body functions independently for two different sublattices A and B? *)
SublatticesQ = False;
(* external AFM field *)
hAFM = 0.0;
(* external pairing field *)
\[CapitalDelta]ext = ConstantArray[0.0, Norb];
(* external staggered on-site field *)
V = 0.0;
(* set True to enforce orbital symmetry and avoid repeating calculations *)
OrbitalSymmetry = False; 
(* enforce P.H. symmetry in the bath *)
ParticleHoleSymmetry = True;
(* list of hoppings for all the orbitals *)
W = ConstantArray[0.25, Norb];

(* add orbital hybridization options if needed ... *)


(*      INPUT PHYSICAL PARAMETERS        *)
(* interaction energy (list of U values for the orbitals) *)
U = {0.1, 0.1}; 
(* Hund's J. It's used only when HundMode = True to enforce rotation invariance of the Kanamori model. *)
JH = 0.00; 
(* density-density opposite spin coupling. It is set automatically if HundMode = True. *)
Ust = 0.1;
(* density-density same spin coupling. It is set automatically if HundMode = True. *)
Usec = 0.1; 
(* pair-hopping coupling. It is set automatically if HundMode = True. *)
Jph = 0.0 * PauliMatrix[1]; 
(* spin-exchange coupling. It is set automatically if HundMode = True. *)
Jse = 0.0; 
(* if this is True, interorbital couplings are authomatically set to Ust=U-2JH; Usec=U-3JH, Jph=JH, Jse=-JH enforcing the rotational invariance *)
HundMode = False; 
(* chemical potential *)
\[Mu] = 0.0; 
(* Crystal field splitting *)
\[Delta] = ConstantArray[0.0, Norb];
(* Raman matrices (one per orbital). Diagonal terms are magnetic field, off-diagonal are Rabi couplings *)
M = {0.1*PauliMatrix[3], 0.0*PauliMatrix[3]};
(* Gauge field associated with Raman tunneling (tensor of dimension Norb x LatticeDim) *)
\[Gamma] = ConstantArray[0.0*Pi*UnitVector[LatticeDim, 1], Norb];
(* temperature *)
T = 0; 


(* INFO ON MATSUBARA AND REAL FREQUENCIES *)
(* fictitious temperature to define Matsubara frequencies *)
TMats = 0.001; 
(* Total number of Matsubara frequencies *)
NMatsubara = 5000; 
(* set min and max value for the set of real frequencies *)
\[Omega]min = -5.0; \[Omega]max = 5.0; 
(* number of real frequencies *)
NReal = 2000;
(* real frequency step *)
d\[Omega] = (\[Omega]max - \[Omega]min)/NReal; 
(* small shift of the pole in the imaginary axis: this avoids singularities, but introduces an artificial broadening of the spectrum *)
\[Eta] = 0.05;


(* INPUT-OUTPUT MANAGEMENT *)
(* directory where the code is *)
CodeDirectory = NotebookDirectory[];
(* directory where the output will be stored *)
OutputDirectory = CodeDirectory<>"prova\\";(*"Jph="<>ToString[Jph[[1,2]]]<>"\\";*)
(* load sectors from a file? *)
LoadSectorsQ = False;
(* load Hamiltonian from a file? *)
LoadHamiltonianQ = False;
(* load symbols from a file *)
LoadSymbolsQ = True;
(* file where relevant sectors are stored *)
SectorsFile = OutputDirectory <> "sectors.m";
(* file where symbols are stored in or imported from *)
SymbolsFile = OutputDirectory <> "symbols.m";


(* NUMERICAL DETAILS OF THE ALGORITHM *)
(* path to input file of bath parameters or "Default" *)
InitializeBathMode = "C:\\Users\\matte\\Desktop\\Mathematica_package\\prova\\hamiltonian_restart.m";
(* if this is True, an extra on-site energy term is added to ensure PH symmetry of the Kanamori Hamiltonian when \[Mu]=0 *)
HFMode = False; 
(* use with caution: if this is True, the full spectrum will be automatically determined! *)
FullDiagonalizationMode = False;
(* below this threshold, two energy levels are assumed to be degenerate *)
DegeneracyThreshold = 10^(-8);
(* dimension threshold above which Lanczos is used to perform exact diagonalization *)
MinLanczosDim = 1000;
(* maximum number of Lanczos iterations for spectrum determination *)
MaxLanczosIter = 2000;
(* if you use Lanczos for spectrum determination, compute AT LEAST this number of eigenstates *)
MinNumberOfEigs = 15;
(* minimum number of Lanczos iterations in the Green function calculation (automatically reduced if it's larger than the Hamiltonian dimension) *)
MinLanczosMomenta = 50;
(* minimum number of DMFT loops *)
DMFTMinIterations = 2; 
(* maximum number of DMFT loops *)
DMFTMaxIterations = 20; 
(* threshold for DMFT loop convergence *)
DMFTerror = 1.0 * 10^(-5); 
(* InverseG0 = Mixing * InverseG0old + (1 - Mixing) * InverseG0 *)
Mixing = 0.5; 
(* type of minimization: "Global" or "Local" *)
MinimizationType = "Local";
(* Method for the minimization procedure *)
MinimizationMethod = "ConjugateGradient";
(* Max number of Conjugate Gradient iterations *)
CGMaxIterations = 2000;
(* Number of Matsubara frequencies used to perform the fit *)
CGNMatsubara = 500;
(* Accuracy goal for minimization procedure *)
CGAccuracy = 6;
(* list of weights attributed to each Matsubara frequency w(i\[Omega]): IT SHOULD BE >= CGNMatsubara *)
CGWeight = ConstantArray[1.0, CGNMatsubara];

(* turn off annoying messages from Eigensystem and stop evaluation after the first error occurs *)
AbortAtFirstErrorQ = False;


Print["Input file imported successfully. "]
