(* ::Package:: *)

(* ::Subtitle:: *)
(*IMPORT DMFT PACKAGE*)


ClearAll["Global`*"];
ClearAll["DMFT`*"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(**)
(*INPUT PARAMETERS FOR DMFT LOOP*)


(*             GENERAL INPUT              *)
Nbath = 2; (* number of bath sites *)
Norb = 1; (* number of orbitals *)
Nimp = 1; (* number of impurity sites *)
L = Nimp + Nbath; (* total number of sites: bath+impurity *)
f = 2; (* number of spin states *)
EdMode = "Normal"; (* call the function EdModeInfo[EdMode] to get details *)

(*      INPUT PHYSICAL PARAMETERS        *)
U = ConstantArray[-0.3, Norb]; (* interaction energy in units of DBethe = 1.0. You have to provide a list of U values in the orbitals *)
JH = 0.0; (* Hund's J. It's used only when HundMode = True to enforce rotation invariance of the Kanamori model. *)
Ust = 0.0; (* density-density opposite spin coupling. It is set automatically if HundMode = True. *)
Usec = 0.0; (* density-density same spin coupling. It is set automatically if HundMode = True. *)
Jph = 0.0; (* pair-hopping coupling. It is set automatically if HundMode = True. *)
Jse = 0.0; (* spin-exchange coupling. It is set automatically if HundMode = True. *)
\[CapitalDelta]0 = 0.2; (* starting value of intraorbital symmetry breaking field if EdMode="Superc" *)
\[CapitalXi]0 = 0.2; (* starting value of interorbital symmetry breaking field if EdMode="InterorbSuperc" or "FullSuperc". *)
\[Mu] = 0; (* chemical potential. It is set automatically if HFMode = True (so you can actually tune it ONLY IF HFMode = False). *)
InitializeBathMode = "Default"; (* path to input file of bath parameters or "Default" *)
HFMode = False; (* if this is True, chemical potential is automatically set to a value that ensures PH symmetry of the Kanamori Hamiltonian *)
HundMode = False; (* if this is True, interorbital couplings are authomatically set to U-2JH and U-3JH, Jph=JH, Jse=-JH enforcing the rotational invariance *)
DegeneracyThreshold = 10^(-9);(* below this threshold, two energy levels are assumed to be degenerate *)

(* OPTIONAL VARIABLES *)
LoadHamiltonianQ = False;(* load Hamiltonian from a file? *)
HnonlocFile = FolderPath<>"Hnonloc_L="<>ToString[L]<>".mx";(* file name for import / export of nonlocal Hamiltonian blocks *)
HlocFile = FolderPath<>"Hloc_L="<>ToString[L]<>".mx";(* file name for import / export of local Hamiltonian blocks *)


(* ::Subtitle:: *)
(**)
(*PREPARATION  (get bath parameters, sectors, Hamiltonians etc...)*)


(* GENERAL VARIABLES *)
LastIteration = False;(*allows to do one more iteration after convergence threshold is reached*)
Converged = False;(*True if DMFT has converged, false otherwise*)
ErrorList = {};(*list of DMFT errors*) 

(* GET BATH AND INTERACTION PARAMETERS *)
Which[
	EdMode == "Normal" || EdMode == "InterorbNormal",
	{e,V} = StartingBath[L, f, Norb, InitializeBathMode, EdMode];
	BathParameters = Flatten@{e,V},
(* ----------------------------------------------------------------- *)
	EdMode == "Superc", 
	{e,V,\[CapitalDelta]} = StartingBath[L, f, Norb, InitializeBathMode, EdMode]; 
	\[CapitalDelta] = \[CapitalDelta]0*\[CapitalDelta]; 
	BathParameters = Flatten@{e, V, \[CapitalDelta]},
(* ----------------------------------------------------------------- *)
	EdMode == "InterorbSuperc",
	{e,V,\[CapitalXi]} = StartingBath[L, f, Norb, InitializeBathMode, EdMode];
	\[CapitalXi] = \[CapitalXi]0*\[CapitalXi];
	BathParameters = Flatten@{e,V,\[CapitalXi]};,
(* ----------------------------------------------------------------- *)
	EdMode == "FullSuperc",
	{e,V,\[CapitalDelta],\[CapitalXi]} = StartingBath[L, f, Norb, InitializeBathMode, EdMode];
	\[CapitalDelta] = \[CapitalDelta]0*\[CapitalDelta]; 
	\[CapitalXi] = \[CapitalXi]0*\[CapitalXi];
	BathParameters = Flatten@{e,V,\[CapitalDelta],\[CapitalXi]};
];
Nparams = Length[BathParameters];(* total number of parameters *)

(* GET INTERACTION PARAMETERS *)
If[
	Norb == 1,
	HundMode = False; Ust = 0.0; Usec = 0.0; Jph = 0.0; Jse = 0.0; JH = 0.0;
];(* avoid stupid bugs: if Norb = 1 set all interorbital couplings to 0 *)
If[
	HundMode,
	Ust = U[[1]] - 2*JH;
	Usec = U[[1]] - 3*JH;
	Jph = JH;
	Jse = -JH;
];(* reset couplings if HundMode = True to enforce rotational invariance *)
If[
	HFMode, 
	\[Mu] = \[Mu] + U[[1]]/2 + (Ust + Usec)*(Norb - 1)/2.;
];(* reset chemical potential if HFMode = True *)
InteractionParameters = Flatten[{U, Ust, Usec, Jph, Jse, \[Mu]}];

(* GET SYMBOLIC GREEN FUNCTION AND WEISS FIELD *)
symbols = Which[
	EdMode == "Normal", 
	Join[
		Table[Symbol["e"<>ToString[i]], {i, L-1}],
		Table[Symbol["V"<>ToString[i]], {i, L-1}]
	],
	EdMode == "Superc",
	Join[
		Table[Symbol["e"<>ToString[i]], {i, L-1}],
		Table[Symbol["V"<>ToString[i]], {i, L-1}],
		Table[Symbol["\[CapitalDelta]"<>ToString[i]], {i, L-1}]
	]
];(* define a suitable list of symbols depending on EdMode *)
G0 = GreenFunction0[L, f, \[Mu], symbols, z, EdMode];
Weiss = Which[
	EdMode == "Normal", FullSimplify[1/G0],
	EdMode == "Superc", FullSimplify[Inverse[G0]]
];
Print["The Weiss field for the given impurity problem is  \!\(\*SuperscriptBox[SubscriptBox[\(G\), \(0\)], \(-1\)]\)(z) = ", Weiss];

(* GET SECTORS *)
QnsSectorList = SectorList[L, f, Norb, EdMode]; (* list of quantum numbers of all the sectors {n,nup} or sz *)
DimSectorList = DimSector[L, f, Norb, #, EdMode]&/@QnsSectorList; (* list of dimensions of all the sectors *)
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList; (* list of all the sectors *)
Print[Style["Recap of input:", 16, Bold]];
Print["Nsectors: ", Length[QnsSectorList], ". Dim. of the largest sector: ", Max@DimSectorList];

(* EXPORT EFFECTIVELY USED INPUT *)
Export[
	FolderPath<>"used_input.dat",{
	{"Nbath =",Nbath}, {"Norb =",Norb}, {"Nimp =",Nimp}, {"f =",f}, {"EdMode =",EdMode}
}];
FilePrint[FolderPath<>"used_input.dat"]

(* GET IMPURITY HAMILTONIAN *)
{HnonlocBlocks, HlocBlocks} = GetHamiltonian[L, f, Norb, Sectors, LoadHamiltonianQ, HnonlocFile, HlocFile, EdMode];


(*                   DMFT LOOP                     *)
(*Do[*)

	(* First print *)
	Print[Style["DMFT Loop n. ", 20, Bold, Red], Style[DMFTiterator, 20, Bold, Red]];
	Print["----------------------------------------------------------------------------------------"];
	Print[Style["        Exact Diagonalization start", 16, Bold, Orange]];
	Print["e = ", e];
	Print["V = ", V];
	If[EdMode=="Superc", Print["\[CapitalDelta] = ", \[CapitalDelta]]];

	(* Build and diagonalize the AIM Hamiltonian + print timing *)
	Print["E.D. time: ", AbsoluteTiming[
	
		Hsectors = HImp[Norb, HnonlocBlocks, HlocBlocks, BathParameters, InteractionParameters, EdMode];
		{EgsSectorList, GsSectorList} = (Map[Eigs[0], Hsectors])\[Transpose];(* find the ground state for each sector *)
		EgsSectorList = Flatten@EgsSectorList;(*correctly reshape the list *)
		GsSectorList = Replace[GsSectorList,{x_List}:>x,{0,-2}];(*correctly reshape the list*)

	]," sec.\n"];
	
	Print["Computing ground state and Green functions..."];

	(* Compute the ground state *)
	Egs = Min@EgsSectorList;(*ground state energy (lowest of all the sectors)*)
	GsSectorIndex =
		Flatten@Position[EgsSectorList,
			_?((Abs[# - Egs] < DegeneracyThreshold)&)
		];(*sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy*)
	DegeneracyWarning = If[Length[GsSectorIndex]>1, True, (*else*) False](*is True if the ground state is degenerate, False otherwise*)
	

(*];*)



