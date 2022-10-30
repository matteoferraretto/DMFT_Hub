(* ::Package:: *)

(* ::Subtitle:: *)
(*IMPORT DMFT PACKAGE*)


ClearAll["Global`*"];
ClearAll["DMFT`*"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"Multiorbital_Package.wl");
<<(FolderPath<>"Facilities.wl");
<<(FolderPath<>"Observables.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*INPUT PARAMETERS FOR DMFT LOOP*)


(*             GENERAL INPUT              *)
Nbath = 4; (* number of bath sites *)
Norb = 2; (* number of orbitals *)
Nimp = 1; (* number of impurity sites *)
L = Nimp + Nbath; (* total number of sites: bath+impurity *)
f = 2; (* number of spin states *)
EdMode = "Normal"; (* call the function EdModeInfo[EdMode] to get details *)
LatticeType = "Bethe"; (* lattice crystal structure: "Bethe", "Hypercubic", etc. *)
LatticeDim = Infinity; (* lattice dimensionality *)
OrbitalSymmetry = True; (* set True to enforce orbital symmetry and avoid repeating calculations *)

(*      INPUT PHYSICAL PARAMETERS        *)
DBethe = ConstantArray[1., Norb]; (* list of half-bandwidths for all the orbitals *)
U = ConstantArray[0.5, Norb]; (* interaction energy in units of DBethe = 1.0. You have to provide a list of U values in the orbitals *)
JH = 0.2; (* Hund's J. It's used only when HundMode = True to enforce rotation invariance of the Kanamori model. *)
Ust = 0.0; (* density-density opposite spin coupling. It is set automatically if HundMode = True. *)
Usec = 0.0; (* density-density same spin coupling. It is set automatically if HundMode = True. *)
Jph = 0.0; (* pair-hopping coupling. It is set automatically if HundMode = True. *)
Jse = 0.0; (* spin-exchange coupling. It is set automatically if HundMode = True. *)
\[Mu] = 0; (* chemical potential. It is set automatically if HFMode = True (so you can actually tune it ONLY IF HFMode = False). *)
T = 0; (* temperature *)
InitializeBathMode = "Default"; (* path to input file of bath parameters or "Default" *)
HFMode = True; (* if this is True, chemical potential is automatically set to a value that ensures PH symmetry of the Kanamori Hamiltonian *)
shift = 0.0; (* energy shift *)
HundMode = True; (* if this is True, interorbital couplings are authomatically set to U-2JH and U-3JH, Jph=JH, Jse=-JH enforcing the rotational invariance *)
DegeneracyThreshold = 10^(-9);(* below this threshold, two energy levels are assumed to be degenerate *)

(* INFO ON MATSUBARA AND REAL FREQUENCIES *)
TMats = 0.001; (* fictitious temperature to define Matsubara frequencies *)
NMatsubara = 5000; (* Total number of Matsubara frequencies *)
i\[Omega] = Table[(2n-1)Pi*I*TMats, {n, NMatsubara}]; (* list of Matsubara frequencies *)
\[Omega]min = -5.; \[Omega]max = 5.; (* set min and max value for the set of real frequencies *)
NReal = 10000;(* number of real frequencies *)
d\[Omega] = (\[Omega]max - \[Omega]min)/NReal; (* real frequency step *)
\[Omega] = Table[\[Omega]min + n*d\[Omega], {n, 0, NReal}]; (* list of real frequencies *)
\[Eta] = 0.025;(* small shift of the pole in the imaginary axis: this avoids singularities, but introduces an artificial broadening of the spectrum *)

(* NUMERICAL DETAILS OF THE ALGORITHM *)
DMFTMinIterations = 2; (* minimum number of DMFT loops *)
DMFTMaxIterations = 20; (* maximum number of DMFT loops *)
DMFTerror = 1.0 * 10^(-5); (* threshold for DMFT loop convergence *)
Mixing = 0.2; (* Mixing * BathParameters + (1 - Mixing) * NewBathParameters *)

(* OPTIONAL VARIABLES *)
LoadHamiltonianQ = True;(* load Hamiltonian from a file? *)
HnonlocFile = FolderPath<>"Hnonloc_L="<>ToString[L]<>"_f="<>ToString[f]<>"_Norb="<>ToString[Norb]<>"_EdMode="<>EdMode<>".mx";(* file name for import / export of nonlocal Hamiltonian blocks *)
HlocFile = FolderPath<>"Hloc_L="<>ToString[L]<>"_f="<>ToString[f]<>"_Norb="<>ToString[Norb]<>"_EdMode="<>EdMode<>".mx";(* file name for import / export of local Hamiltonian blocks *)


(* ::Subtitle:: *)
(*PREPARATION  (get bath parameters, sectors, Hamiltonians etc...)*)


(* GENERAL VARIABLES *)
LastIteration = False; (* allows to do one more iteration after convergence threshold is reached *)
Converged = False; (* True if DMFT has converged, false otherwise *)
ErrorList = {}; (* list of DMFT errors *) 

(* GET BATH PARAMETERS *)
BathParameters = StartingBath[L, f, Norb, InitializeBathMode, EdMode];
Nparams = Length[BathParameters];(* total number of parameters *)

(* GET INTERACTION PARAMETERS *)
If[Norb == 1,
	HundMode = False; Ust = 0.0; Usec = 0.0; Jph = 0.0; Jse = 0.0; JH = 0.0;
];(* avoid stupid bugs: if Norb = 1 set all interorbital couplings to 0 *)
If[HundMode,
	Ust = U[[1]] - 2*JH; Usec = U[[1]] - 3*JH; Jph = 0*JH; Jse = -0*JH;
];(* reset couplings if HundMode = True to enforce rotational invariance *)
If[HFMode, 
	shift = -U[[1]]/2 - (Ust + Usec)*(Norb - 1)/2.;
];(* reset chemical potential if HFMode = True *)
InteractionParameters = Flatten[{U, Ust, Usec, Jph, Jse, \[Mu] + shift}];

(* GET SYMBOLIC GREEN FUNCTION AND WEISS FIELD *)
symbols = Symbols[L, f, EdMode]; (* define a suitable list of symbols depending on EdMode *)
Weiss = Apart[WeissField[L, f, \[Mu], symbols, z, EdMode], z];
Print["The Weiss field for the given impurity problem is  \!\(\*SuperscriptBox[SubscriptBox[\(G\), \(0\)], \(-1\)]\)(z) = ", Weiss];

(* GET SECTORS *)
QnsSectorList = SectorList[L, f, Norb, EdMode]; (* list of quantum numbers of all the sectors {n,nup} or sz *)
DimSectorList = DimSector[L, f, Norb, #, EdMode]&/@QnsSectorList; (* list of dimensions of all the sectors *)
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList; (* list of all the sectors *)
SectorsDispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&, QnsSectorList], 1]]; (* list of rules that assigns every element of QnsSectorList to an integer *)
Print[Style["Recap of input:", 16, Bold]];
Print["Nsectors: ", Length[QnsSectorList], ". Dim. of the largest sector: ", Max@DimSectorList];

(* EXPORT EFFECTIVELY USED INPUT *)
Export[
	FolderPath<>"used_input.dat",{
	{"Nbath =",Nbath}, {"Norb =",Norb}, {"Nimp =",Nimp}, {"f =",f}, {"EdMode =",EdMode}, {"U =",U}
}];
FilePrint[FolderPath<>"used_input.dat"]

(* GET IMPURITY HAMILTONIAN *)
{HnonlocBlocks, HlocBlocks} = GetHamiltonian[L, f, Norb, Sectors, LoadHamiltonianQ, HnonlocFile, HlocFile, EdMode];


(*                   DMFT LOOP                     *)
Do[
	ClearAll[Hsectors, EgsSectorList, GsSectorList];
	
	(* First print *)
	Print[Style["DMFT Loop n. ", 20, Bold, Red], Style[DMFTiterator, 20, Bold, Red]];
	Print[Style["\t\t Exact Diagonalization start", 16, Bold, Orange]];
	Print["Bath parameters: ", BathParameters];

	(* Build and diagonalize the AIM Hamiltonian + print timing *)
	Print["E.D. time: ", AbsoluteTiming[
	
		Hsectors = HImp[Norb, HnonlocBlocks, HlocBlocks, BathParameters, InteractionParameters, EdMode];
		eigs = Map[Eigs[#, "Temperature" -> T, "MinLanczosDim"->100, "MaxIterations" -> 2000]&, Hsectors];
		EgsSectorList = eigs[[All, 1]];
		GsSectorList = eigs[[All, 2]];
		ClearAll[eigs];
		
	]," sec.\n"];
	
	
	(* ZERO TEMPERATURE CALCULATIONS *)
	If[T == 0,
		Print["Computing ground state and Green functions..."];
		(* Compute the ground state *)
		Egs = Min[Flatten[EgsSectorList]];(* ground state energy (lowest of all the sectors) *)
		GsSectorIndex = Position[
			EgsSectorList,
			_?((Abs[# - Egs] < DegeneracyThreshold)&)
		];(* sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy *)
		Gs = MapApply[GsSectorList[[##]]&, GsSectorIndex];(* list of all the degenerate ground states *)
		GsQns = QnsSectorList[[GsSectorIndex[[All, 1]]]];(* list of quantum numbers of the degenerate ground states *)	
		Print["\t\t Ground state info:\n", "Egs = ", Egs, "   Quantum numbers = ", GsQns];(* print relevant information about the ground state *)
		
		Do[
			Print["impurity density [orb="<>ToString[orb]<>"] = ", Sum[Density[L, f, Norb, 1, \[Sigma], orb, Sectors, EgsSectorList, GsSectorList, T], {\[Sigma], f}] ];
			Print["impurity double occupancy [orb="<>ToString[orb]<>"] = ", SquareDensity[L, f, Norb, {1,1}, {1,2}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T] ];
			If[EdMode=="Superc", Print["order parameter [orb="<>ToString[orb]<>"] = ", CdgCdg[L, f, Norb, {1,1}, {1,2}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T] ];]
		, {orb, Norb}];
		
		If[OrbitalSymmetry,
			IndependentParameters = TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode];
			(* G^-1(i\[Omega]) *)
			InverseG = Mean[MapApply[
				InverseGreenFunction[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
				{Gs, GsQns}\[Transpose]
			]];
			(* Subscript[G, 0]^-1(i\[Omega]) *)
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0];
			InverseG0 = (Weiss/.Thread[symbols -> IndependentParameters])/.{z -> #}&/@i\[Omega];
			(* \[CapitalSigma](i\[Omega]) *)
			\[CapitalSigma] = InverseG0 - InverseG;
			(* Subscript[G, loc](i\[Omega]) *)
			LocalG = LocalGreenFunction[DBethe[[1]], \[Mu], \[CapitalSigma], EdMode, i\[Omega], Lattice -> LatticeType, LatticeDimension -> LatticeDim, NumberOfPoints -> 3000];
			(* Self consistency *)
			Print["Quasiparticle weight z = ", QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode] ];
			Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
			Print["S.C. time: ", First@AbsoluteTiming[
			
			NewBathParameters = ReshapeBathParameters[L, f, Norb,	
				SelfConsistency[DBethe[[1]], \[Mu], Weiss, symbols, z, IndependentParameters, LocalG, \[CapitalSigma], i\[Omega], EdMode, 
				Lattice -> LatticeType, LatticeDimension -> LatticeDim, Minimum -> "Local", NumberOfFrequencies -> 500, MaxIterations -> 2000, AccuracyGoal -> 7],
			OrbitalSymmetry, EdMode];
			
			], " sec." ];
			error = DMFTError[InverseG0[[;;2000]], InverseG0old[[;;2000]], EdMode];,
		(* --------------------------  *)	
		(* else, if no orbital symmetry: *)
			IndependentParameters = Table[
				TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode],
			{orb, Norb}];
			(* { G^-1Subscript[(i\[Omega]), orb=1] , G^-1Subscript[(i\[Omega]), orb=2] , ...} *)
			InverseG = Table[
				Mean[MapApply[
					InverseGreenFunction[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
					{Gs, GsQns}\[Transpose]
				]], {orb, Norb}];
			(* { Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=1] , Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=2] , ...} *)
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0];
			InverseG0 = Table[
				(Weiss/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode]])/.{z -> #}&/@i\[Omega],
				{orb, Norb}];
			(* { \[CapitalSigma]Subscript[(i\[Omega]), orb=1] , \[CapitalSigma]Subscript[(i\[Omega]), orb=2] , ...} *)
			\[CapitalSigma] = InverseG0 - InverseG;
			(* { Subscript[G, loc]Subscript[(i\[Omega]), orb=1] , Subscript[G, loc]Subscript[(i\[Omega]), orb=2] , ...} *)
			LocalG = Table[
				LocalGreenFunction[DBethe[[orb]], \[Mu], \[CapitalSigma][[orb]], EdMode, i\[Omega], Lattice -> LatticeType, LatticeDimension -> LatticeDim, NumberOfPoints -> 3000]
			, {orb, Norb}];
			(* Self consistency *)
			Print["Quasiparticle weight z = ", Table[QuasiparticleWeight[\[CapitalSigma][[orb]], i\[Omega], EdMode], {orb, Norb}] ];
			Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
			Print["S.C. time: ", First@AbsoluteTiming[
			
			NewBathParameters = ReshapeBathParameters[L, f, Norb, Table[
				SelfConsistency[DBethe[[orb]], \[Mu], Weiss, symbols, z, IndependentParameters[[orb]], LocalG[[orb]], \[CapitalSigma][[orb]], i\[Omega], EdMode,
				Lattice -> LatticeType, LatticeDimension -> LatticeDim, Minimum -> "Local", NumberOfFrequencies -> 500, MaxIterations -> 2000, AccuracyGoal -> 7]
			, {orb, Norb}], OrbitalSymmetry, EdMode];
			
			], " sec." ];
			error = (1./Norb)*Sum[
				DMFTError[InverseG0[[orb]], InverseG0old[[orb]], EdMode],
			{orb, Norb}];
		];
	];
	
	(* FINITE TEMPERATURE CALCULATIONS *)
	If[T != 0,
		Print["not supported."];
		Break[];
	];
	
	(* update bath parameters *)
	BathParameters = Mixing * BathParameters + (1 - Mixing) * NewBathParameters;
	Print["DMFT error: ", error];
	
	(*Print*)
	Print[Style["\t\t Self Consistency completed", 16, Bold, Magenta]];
	Print["----------------------------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------------------------"];

	(* Exit DMFT Loop if convengerce is reached *)
	If[DMFTiterator > DMFTMinIterations && error < DMFTerror && LastIteration, Break[];];
	If[error < DMFTerror, LastIteration = True, (*else*) LastIteration = False];


, {DMFTiterator, DMFTMaxIterations}]



PlotMatsubara[Im[\[CapitalSigma]], i\[Omega], EdMode]
PlotMatsubara[Re[\[CapitalSigma]], i\[Omega], EdMode]
(*
ListPlot[{Abs[i\[Omega]],Re[\[CapitalSigma][[All,1,2]]]}\[Transpose],Joined->True, PlotRange->{0,4}]
*)
spectralfunction = SpectralFunction[L, f, Norb, 1, 1, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega], \[Eta]];
ListPlot[spectralfunction, Joined->True, PlotRange->All]
d\[Omega] * Total[spectralfunction[[All,2]]]
(*
With[
	{G = Inverse[#] &/@ InverseG},
	- TMats * Total[G[[All, 1, 2]]]
]*)


ListPlot[{
	Im[1./LocalG],
	Im[InverseG]
	}, Joined->True, PlotRange->{{0,5000},Automatic}, PlotStyle->{Thick, Dashing[.05]}]

ListPlot[Im[\[CapitalSigma]], Joined->True]

\[Eta] = 0.005
SelfEnergyRealFreq = (Weiss/.Thread[symbols -> IndependentParameters])/.{z -> #}&/@(\[Omega]+I*\[Eta]) - 1./Mean[MapApply[
	GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega]+I*\[Eta]]&,
	{Gs, GsQns}\[Transpose]
]]

ListPlot[{\[Omega], Re[SelfEnergyRealFreq]}\[Transpose], Joined->True]


\[Epsilon][kx_, ky_] := -2(Cos[kx] + Cos[ky]);
(*pos[\[Omega]_] := Floor[(\[Omega] - \[Omega]min)/d\[Omega]] + 1;*)
spectralWeight[kx_, ky_, pos_] := 1./(\[Omega][[pos]] + I*\[Eta] - \[Epsilon][kx,ky] - SelfEnergyRealFreq[[pos]])

spectraldata = Table[{k, \[Omega][[pos]], -(1./Pi)*Im[spectralWeight[k, 0, pos]]}, {k,-Pi,Pi,0.05Pi}, {pos, Range[1, NReal-4000]}];
Dimensions[spectraldata]


Show[
ListDensityPlot[Transpose@Partition[Flatten[spectraldata,1], 21]],
Plot[\[Epsilon][k,0], {k,-Pi,Pi}]
]
