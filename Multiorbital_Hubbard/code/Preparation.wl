(* ::Package:: *)

(* GENERAL VARIABLES *)
(* allows to do one more iteration after convergence threshold is reached *)
LastIteration = False; 
(* True if DMFT has converged, false otherwise *)
Converged = False; 
(* initialize list of DMFT errors *) 
ErrorList = {}; 


(* GET BATH PARAMETERS *)
BathParameters = StartingBath[L, f, Norb, InitializeBathMode, EdMode];
Nparams = Length[BathParameters];


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
	shift = -U[[1]]/2 - (Ust + Usec)*(Norb - 1)/2.;
];
InteractionParameters = Flatten[{U, Ust, Usec, Jph, Jse, \[Mu] + shift}];


(* GET SYMBOLIC GREEN FUNCTION AND WEISS FIELD *)
(* define a suitable list of symbols depending on EdMode *)
symbols = Symbols[L, f, Norb, EdMode]; 
(* define a symbolic expression for the Weiss field *)
Weiss = WeissField[L, f, Norb, \[Mu], symbols, z, EdMode];
(* print it *)
Print["The Weiss field for the given impurity problem is  \!\(\*SuperscriptBox[SubscriptBox[\(G\), \(0\)], \(-1\)]\)(z) = ", Weiss];


(* GET SECTORS *)
(* list of quantum numbers of all the sectors *)
QnsSectorList = SectorList[L, f, Norb, EdMode]; 
(* list of dimensions of all the sectors *)
DimSectorList = DimSector[L, f, Norb, #, EdMode]&/@QnsSectorList; 
(* list of all the sectors *)
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList; 
(* list of rules that assigns every element of QnsSectorList to an integer *)
SectorsDispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&, QnsSectorList], 1]]; 
(* print recap *)
Print[Style["Recap of input:", 16, Bold]];
Print["Nsectors: ", Length[QnsSectorList], ". Dim. of the largest sector: ", Max@DimSectorList];


(* EXPORT EFFECTIVELY USED INPUT *)
(* Export[
	OutputDirectory<>"used_input.dat",{
	{"Nbath =",Nbath}, {"Norb =",Norb}, {"Nimp =",Nimp}, {"f =",f}, {"EdMode =",EdMode}, {"U =",U}
}];
FilePrint[OutputDirectory<>"used_input.dat"] *)


(* GET LATTICE *)
{LatticeEnergies, LatticeWeights} = GetLatticeEnergies[W, LatticeType, LatticeDim, LatticePoints];


(* GET IMPURITY HAMILTONIAN *)
{HnonlocBlocks, HlocBlocks} = GetHamiltonian[L, f, Norb, Sectors, LoadHamiltonianQ, HnonlocFile, HlocFile, EdMode];
