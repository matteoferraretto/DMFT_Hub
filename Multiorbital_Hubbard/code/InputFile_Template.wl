(* ::Package:: *)

(*             GENERAL INPUT              *)
(* number of bath sites *)
Nbath = 2; 
(* number of orbitals *)
Norb = 2; 
(* number of impurity sites *)
Nimp = 1; 
(* total number of sites: bath+impurity *)
L = Nimp + Nbath; 
(* number of spin states *)
f = 2; 
(* call the function EdModeInfo[EdMode] to get details *)
EdMode = "FullSuperc"; 
(* lattice crystal structure: "Bethe", "Hypercubic", etc. *)
LatticeType = "Hypercubic"; 
(* lattice dimensionality *)
LatticeDim = 2; 
(* lattice points *)
LatticePoints = 1600;
(* set True to enforce orbital symmetry and avoid repeating calculations *)
OrbitalSymmetry = True; 
(* list of half-bandwidths for all the orbitals *)
W = ConstantArray[1., Norb]; 
(* add orbital hybridization options if needed ... *)


(*      INPUT PHYSICAL PARAMETERS        *)
(* interaction energy (list of U values for the orbitals) *)
U = ConstantArray[-0.3, Norb]; 
(* Hund's J. It's used only when HundMode = True to enforce rotation invariance of the Kanamori model. *)
JH = 0.0; 
(* density-density opposite spin coupling. It is set automatically if HundMode = True. *)
Ust = 0.0; 
(* density-density same spin coupling. It is set automatically if HundMode = True. *)
Usec = 0.0; 
(* pair-hopping coupling. It is set automatically if HundMode = True. *)
Jph = 0.0; 
(* spin-exchange coupling. It is set automatically if HundMode = True. *)
Jse = 0.0; 
(* if this is True, interorbital couplings are authomatically set to Ust=U-2JH; Usec=U-3JH, Jph=JH, Jse=-JH enforcing the rotational invariance *)
HundMode = True; 
(* chemical potential *)
\[Mu] = 0; 
(* temperature *)
T = 0; 


(* INFO ON MATSUBARA AND REAL FREQUENCIES *)
(* fictitious temperature to define Matsubara frequencies *)
TMats = 0.001; 
(* Total number of Matsubara frequencies *)
NMatsubara = 5000; 
(* list of Matsubara frequencies (don't touch) *)
i\[Omega] = Table[(2n-1)Pi*I*TMats, {n, NMatsubara}]; 
(* set min and max value for the set of real frequencies *)
\[Omega]min = -5.; \[Omega]max = 5.; 
(* number of real frequencies *)
NReal = 10000;
(* real frequency step *)
d\[Omega] = (\[Omega]max - \[Omega]min)/NReal; 
(* list of real frequencies *)
\[Omega] = Table[\[Omega]min + n*d\[Omega], {n, 0, NReal}]; 
(* small shift of the pole in the imaginary axis: this avoids singularities, but introduces an artificial broadening of the spectrum *)
\[Eta] = 0.025;


(* NUMERICAL DETAILS OF THE ALGORITHM *)
(* path to input file of bath parameters or "Default" *)
InitializeBathMode = "Default"; 
(* if this is True, an extra on-site energy term is added to ensure PH symmetry of the Kanamori Hamiltonian when \[Mu]=0 *)
HFMode = True; 
(* energy shift *)
shift = 0.0; 
(* below this threshold, two energy levels are assumed to be degenerate *)
DegeneracyThreshold = 10^(-9);
(* minimum number of DMFT loops *)
DMFTMinIterations = 2; 
(* maximum number of DMFT loops *)
DMFTMaxIterations = 20; 
(* threshold for DMFT loop convergence *)
DMFTerror = 1.0 * 10^(-5); 
(* Mixing * BathParameters + (1 - Mixing) * NewBathParameters *)
Mixing = 0; 


(* INPUT-OUTPUT MANAGEMENT *)
(* directory where the code is *)
CodeDirectory = NotebookDirectory[];
(* directory where the output will be stored *)
OutputDirectory = CodeDirectory<>"U="<>ToString[U[[1]]]<>"\\";
(* load Hamiltonian from a file? *)
LoadHamiltonianQ = False;
(* file name for import / export of nonlocal Hamiltonian blocks *)
HnonlocFile = OutputDirectory<>"Hnonloc_L="<>ToString[L]<>"_f="<>ToString[f]<>"_Norb="<>ToString[Norb]<>"_EdMode="<>EdMode<>".mx";
(* file name for import / export of local Hamiltonian blocks *)
HlocFile = OutputDirectory<>"Hloc_L="<>ToString[L]<>"_f="<>ToString[f]<>"_Norb="<>ToString[Norb]<>"_EdMode="<>EdMode<>".mx";



Print["Input file imported successfully. "]
