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
Nbath = 4; (*number of bath sites*)
Norb = 2; (*number of orbitals*)
Nimp = 1; (*number of impurity sites*)
L = Nimp + Nbath; (*total number of sites: bath+impurity*)
f = 2; (* number of spin states *)
EdMode = "InterorbSuperc"; (* "Normal" = no symmetry breaking;  "Superc" = bath exchanging pairs with reservoir *)

(*      INPUT PHYSICAL PARAMETERS        *)
U = -0.3; (* interaction energy in units of DBethe = 1.0 *)
InitializeBathMode = "Default"; (*path to input file of bath parameters or "Default"*)
\[CapitalDelta]0 = 0.2; (* starting value of symmetry breaking field if EdMode="Superc" *)
\[CapitalXi]0 = 0.2; (* starting value of interorbital symmetry breaking field if EdMode="Superc" *)
\[Mu] = 0; (*chemical potential*)


(* ::Subtitle:: *)
(**)
(*PREPARATION  (get bath parameters, sectors, Hamiltonians etc...)*)


(* GENERAL VARIABLES *)
LastIteration = False;(*allows to do one more iteration after convergence threshold is reached*)
Converged = False;(*True if DMFT has converged, false otherwise*)
ErrorList = {};(*list of DMFT errors*) 

(* GET BATH PARAMETERS *)
Which[
	EdMode == "Normal" || EdMode == "InterorbNormal",
	{e,V} = StartingBath[InitializeBathMode, Nbath, Norb, EdMode];
	Parameters = Flatten@{e,V},
(* ----------------------------------------------------------------- *)
	EdMode == "Superc", 
	{e,V,\[CapitalDelta]} = StartingBath[InitializeBathMode, Nbath, Norb, EdMode]; 
	\[CapitalDelta] = \[CapitalDelta]0*\[CapitalDelta]; 
	Parameters = Flatten@{e, V, \[CapitalDelta]},
(* ----------------------------------------------------------------- *)
	EdMode == "InterorbSuperc",
	{e,V,\[CapitalDelta],\[CapitalXi]} = StartingBath[InitializeBathMode, Nbath, Norb, EdMode];
	\[CapitalDelta] = \[CapitalDelta]0*\[CapitalDelta]; 
	\[CapitalXi] = \[CapitalXi]0*\[CapitalXi];
	Parameters = Flatten@{e,V,\[CapitalDelta],\[CapitalXi]};
];
Nparams = Length[Parameters];(* total number of parameters *)

(* GET SECTORS *)
QnsSectorList = SectorList[L, f, Norb, EdMode]; (* list of quantum numbers of all the sectors {n,nup} or sz *)
DimSectorList = DimSector[L, f, Norb, #, EdMode]&/@QnsSectorList; (* list of dimensions of all the sectors *)
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList; (* list of all the sectors *)
SectorsDispatch = Dispatch@Flatten[MapIndexed[{#1->#2[[1]]}&, QnsSectorList], 1]; (* dispatch of the sectors, i.e. rule that associates an integer to some quantum numbers *)

(* SECTORS RECAP *)
Print[Style["Recap of input:", 16, Bold]];
Print["Nsectors: ", Length[QnsSectorList], ". Dim. of the largest sector: ", Max@DimSectorList];



