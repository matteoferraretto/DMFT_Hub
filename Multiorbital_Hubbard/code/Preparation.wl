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


(* AVOID BUGS AND SHOW WARNINGS *)
(* set OrbitalSymmetry in some cases *)
If[EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc", OrbitalSymmetry = False];
If[Norb == 1, OrbitalSymmetry = True];
(* avoid stupid bugs related to the C.G. *)
If[CGNMatsubara > NMatsubara, CGNMatsubara = NMatsubara];
If[Length[CGWeight] != CGNMatsubara, 
	Print[Style["Error. CGWeight should be a list of length " <> ToString[CGNMatsubara] <> ", proceeding with default options.", Red]];
	CGWeight = ConstantArray[1., CGNMatsubara];
];
(* avoid stupid bugs related to gauge field in Raman *)
If[M == ConstantArray[0.0, {Norb, f, f}] && EdMode == "Raman",
	Print[Style["Warning. The Raman matrix is zero but EdMode = ''Raman''. You should either add a small Raman coupling or change EdMode to ''Normal''.", Red]];
];
Do[
	If[Length[M[[orb]]] != f && EdMode == "Raman", 
		Print[Style["Error. The Raman matrix associated to orbital "<>ToString[orb]<>" has a wrong dimension.", Red]];
		Abort[];
	];
	If[Norm[\[Gamma][[orb]]] != 0.0 && LatticeType == "Bethe" && EdMode == "Raman",
		Print[Style["Error. Gauge flux is not compatible with Bethe lattice due to lack of geometrical structue. Proceeding with \[Gamma] = 0.", Red]];
		\[Gamma][[orb]] = 0.0;
	];
	If[Norm[\[Gamma][[orb]]] != 0.0 && Length[\[Gamma][[orb]]] != LatticeDim && EdMode == "Raman",
		Print[Style["Error. The vector \[Gamma] should have the same dimension as the lattice.", Red]];
		Abort[];
	];
	If[M[[orb]] == 0.0*M[[orb]] && \[Gamma][[orb]] != 0.0*\[Gamma][[orb]],
		Print[Style["Warning. In orbital "<>ToString[orb]<>" there's no Raman field, but a gauge field has been provided. ", Red]];
	];
, {orb, Norb}];
If[EdMode != "Raman" && M != ConstantArray[0.0, {Norb,f,f}],
	Print[Style["Error. A Raman matrix is provided, but EdMode is not set to ''Raman''. Switch EdMode to ''Raman'' or remove the Raman field.", Red]];
	Abort[];
];
(* avoid bugs related to sublattice calculation mode *)
If[!SublatticesQ && (V != 0.0||hAFM != 0.0),
	Print[Style["Warning. The staggered external field is different from 0.0, but SublatticeQ = False. Proceeding setting it to 0.", Red]];
	V = 0.0; hAFM = 0.0;
];
(* sublattice calculation only supported for 2 spin states *)
If[SublatticesQ && EdMode == "Raman" && f>2,
	Print[Style["Error. Antiferromagnetic calculation is not supported for more than 2 flavors.", Red]];
	Abort[];
];


(* GET INTERACTION PARAMETERS *)
(* avoid stupid bugs: if Norb = 1 set all interorbital couplings to 0 *)
If[Norb == 1,
	HundMode = False; Ust = 0.0; Usec = 0.0; Jse = 0.0; JH = 0.0;
];
(* reset couplings if HundMode = True to enforce rotational invariance *)
If[HundMode,
	Ust = U[[1]] - 2*JH; 
	Usec = U[[1]] - 3*JH; 
	Jph = 0 * JH * (ConstantArray[1, {Norb, Norb}] - IdentityMatrix[Norb]); (* just off diagonal pair hopping *)
	Jse = - 0 * JH;
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


(* GET BATH AND INTERACTION PARAMETERS *)
BathParameters = StartingBath[L, f, Norb, \[Delta]-\[Mu], InitializeBathMode, EdMode, SublatticesQ, V0 -> 1.0, \[CapitalDelta]0 -> 0.2, \[CapitalXi]0 -> 0.2, \[CapitalOmega]0 -> 0.2];
InteractionParameters = LocalParameters[M, \[Delta], U, Ust, Usec, Jph, Jse, \[Mu], shift, SublatticesQ, hAFM, V];
symbols = Symbols[L, f, Norb, EdMode]; (* define a suitable list of symbols depending on EdMode *)


(* GET SECTORS *)
(* list of quantum numbers of all the sectors *)
QnsSectorList = SectorList[L, f, Norb, EdMode]; 
(* list of dimensions of all the sectors *)
DimSectorList = DimSector[L, f, Norb, #, EdMode] &/@ QnsSectorList; 
(* list of all the sectors *)
Sectors = BuildSector[L, f, Norb, #, EdMode] &/@ QnsSectorList; 
(* list of rules that assigns every element of QnsSectorList to an integer *)
SectorsDispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&, QnsSectorList], 1]]; 
(* if full diagonalization is required, make sure that Lanczos will never be used *)
If[FullDiagonalizationMode, MinLanczosDim = Max @ DimSectorList];
(* print recap *)
Print[Style["Recap of input:", 16, Bold]];
Print["Nsectors: ", Length[QnsSectorList], ". Dim. of the largest sector: ", Max @ DimSectorList];


(* GET LATTICE ENERGIES *)
If[EdMode == "Raman",
	If[SublatticesQ,
		{LatticeEnergies, LatticeWeights} = GetLatticeEnergiesRamanSublattices[W, \[Delta], M, \[Gamma], LatticeType, LatticeDim, LatticePoints],
		(* else *)
		{LatticeEnergies, LatticeWeights} = GetLatticeEnergiesRaman[W, \[Delta], M, \[Gamma], LatticeType, LatticeDim, LatticePoints]
	],
	(* else *)
	If[SublatticesQ,
		{LatticeEnergies, LatticeWeights} = GetLatticeEnergiesSublattices[W, \[Delta], LatticeType, LatticeDim, LatticePoints],
		(* else *)
		{LatticeEnergies, LatticeWeights} = GetLatticeEnergies[W, \[Delta], LatticeType, LatticeDim, LatticePoints]
	];
];


(* GET IMPURITY HAMILTONIAN *)
(* file name for import / export of nonlocal Hamiltonian blocks *)
HnonlocFile = CodeDirectory<>"Hnonloc_L="<>ToString[L]<>"_f="<>ToString[f]<>"_Norb="<>ToString[Norb]<>"_EdMode="<>EdMode<>".mx";
(* file name for import / export of local Hamiltonian blocks *)
HlocFile = CodeDirectory<>"Hloc_L="<>ToString[L]<>"_f="<>ToString[f]<>"_Norb="<>ToString[Norb]<>"_EdMode="<>EdMode<>".mx";
(* get hamiltonians *)
{HnonlocBlocks, HlocBlocks} = GetHamiltonian[L, f, Norb, Nimp, Sectors, LoadHamiltonianQ, HnonlocFile, HlocFile, EdMode];


(* COPY INPUT FILE *)
Save[
	OutputDirectory<>"InputFile_Used.wl",
	{Nbath, Norb, Nimp, L, f, EdMode, LatticeType, LatticeDim, LatticePoints, SublatticesQ, V, OrbitalSymmetry, W, U, JH, Ust, Usec, Jph, Jse, HundMode, \[Mu], \[Delta], T, TMats, NMatsubara, 
	\[Omega]min, \[Omega]max, NReal, d\[Omega], \[Eta], CodeDirectory, OutputDirectory, LoadHamiltonianQ, InitializeBathMode, HFMode, FullDiagonalizationMode, DegeneracyThreshold, MinLanczosDim, MaxLanczosIter, 
	MinNumberOfEigs, DMFTMinIterations, DMFTMaxIterations, DMFTerror, Mixing, MinimizationType, MinimizationMethod, CGMaxIterations, CGNMatsubara, CGAccuracy, CGWeight}
]


(* TURN OFF ANNOYING MESSAGES AND ABORT AS SOON AS AN ERROR SHOWS UP *)
If[AbortAtFirstErrorQ,
	Off[Eigensystem::arh];
	messageHandler = If[Last[#], Abort[]] &;
	Internal`AddHandler["Message", messageHandler]
]
