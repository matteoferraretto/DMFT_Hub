(* ::Package:: *)

(* GENERAL VARIABLES *)
(* allows to do one more iteration after convergence threshold is reached *)
LastIteration = False; 
(* True if DMFT has converged, false otherwise *)
Converged = False; 
(* initialize list of DMFT errors *) 
ErrorList = {}; 
(* initialize Matsubara frequencies [use just CGNMatsubara frequencies, not NMatsubara!] *)
i\[Omega] = Table[(2n-1)Pi*I*TMats, {n, CGNMatsubara}]; 

(* set OrbitalSymmetry to False in some cases to avoid stupid bugs *)
If[EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc", OrbitalSymmetry = False];
(* avoid a bunch of stupid bugs related to the C.G. *)
If[CGNMatsubara > NMatsubara, CGNMatsubara = NMatsubara];
If[Length[CGWeight] != CGNMatsubara, 
	Print[Style["Error. CGWeight should be a list of length " <> ToString[CGNMatsubara] <> ", proceeding with default options.", Red]];
	CGWeight = ConstantArray[1., CGNMatsubara];
];


(* GET INTERACTION PARAMETERS *)
(* avoid stupid bugs: if Norb = 1 set all interorbital couplings to 0 *)
If[Norb == 1,
	HundMode = False; Ust = 0.0; Usec = 0.0; Jph = 0.0; Jse = 0.0; JH = 0.0;
];
(* reset couplings if HundMode = True to enforce rotational invariance *)
If[HundMode,
	Ust = U[[1]] - 2*JH; Usec = U[[1]] - 3*JH; Jph = 0*JH; Jse = -0*JH;
];
(* shift on site energy if HFMode = True *)
If[HFMode, 
	shift = -U[[1]]/2 - (Ust + Usec)*(Norb - 1)/2.;,
(* else *)
	shift = 0;
];
(* avoid stupid bugs: check if \[Delta] has the correct length and that there is no conflict with Orbital symmetry requirement *)
If[Length[\[Delta]] != Norb, 
	Print[Style["Error. Crystal field splitting should be a list of " <> ToString[Norb] <> " elements. Proceeding with no crystal field splitting.", Red]];
	\[Delta] = ConstantArray[0.0, Norb];
];
If[OrbitalSymmetry && Norb > 1,
	If[\[Delta] != ConstantArray[0.0, Norb],
		Print[Style["Error. Orbital symmetry is incompatible with a crystal field splitting. Proceeding with no crystal field splitting.", Red]];
	];
	\[Delta] = ConstantArray[0.0, Norb];
];

(* get flat list of interaction parameters *)
InteractionParameters = Flatten[{\[Delta], U, Ust, Usec, Jph, Jse, - \[Mu] + shift}];


(* GET BATH PARAMETERS *)
BathParameters = StartingBath[L, f, Norb, \[Delta]-\[Mu], InitializeBathMode, EdMode, V0 -> 0.1, \[CapitalDelta]0 -> 1.0, \[CapitalXi]0 -> 1.0];
Nparams = Length[BathParameters];


(* GET SYMBOLIC GREEN FUNCTION AND WEISS FIELD *)
(* define a suitable list of symbols depending on EdMode *)
symbols = Symbols[L, f, Norb, EdMode]; 
(* define a symbolic expression for the Weiss field *)
(* Weiss = WeissField[L, f, Norb, \[Mu]eff, symbols, z, EdMode]; 
(* print it *)
If[EdMode != "FullSuperc",
	Print["The Weiss field for the given impurity problem is  \!\(\*SuperscriptBox[SubscriptBox[\(G\), \(0\)], \(-1\)]\)(z) = ", Weiss];
];*)


(* GET SECTORS *)
(* list of quantum numbers of all the sectors *)
QnsSectorList = SectorList[L, f, Norb, EdMode]; 
(* list of dimensions of all the sectors *)
DimSectorList = DimSector[L, f, Norb, #, EdMode]&/@QnsSectorList; 
(* list of all the sectors *)
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList; 
(* list of rules that assigns every element of QnsSectorList to an integer *)
SectorsDispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&, QnsSectorList], 1]]; 
(* if full diagonalization is required, make sure that Lanczos will never be used *)
If[FullDiagonalizationMode, MinLanczosDim = Max @ DimSectorList];
(* print recap *)
Print[Style["Recap of input:", 16, Bold]];
Print["Nsectors: ", Length[QnsSectorList], ". Dim. of the largest sector: ", Max @ DimSectorList];


(* GET LATTICE *)
{LatticeEnergies, LatticeWeights} = GetLatticeEnergies[W, \[Delta], LatticeType, LatticeDim, LatticePoints];


(* GET IMPURITY HAMILTONIAN *)
(* file name for import / export of nonlocal Hamiltonian blocks *)
HnonlocFile = OutputDirectory<>"Hnonloc_L="<>ToString[L]<>"_f="<>ToString[f]<>"_Norb="<>ToString[Norb]<>"_EdMode="<>EdMode<>".mx";
(* file name for import / export of local Hamiltonian blocks *)
HlocFile = OutputDirectory<>"Hloc_L="<>ToString[L]<>"_f="<>ToString[f]<>"_Norb="<>ToString[Norb]<>"_EdMode="<>EdMode<>".mx";
(* get hamiltonians *)
{HnonlocBlocks, HlocBlocks} = GetHamiltonian[L, f, Norb, Nimp, Sectors, LoadHamiltonianQ, HnonlocFile, HlocFile, EdMode];


(* COPY INPUT FILE *)
CopyFile[
	CodeDirectory<>"InputFile_Template.wl",
	OutputDirectory<>"InputFile_Used.wl",
	OverwriteTarget -> True
];
