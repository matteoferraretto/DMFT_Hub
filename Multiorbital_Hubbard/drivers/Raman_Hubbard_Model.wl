(* ::Package:: *)

$Path = Join[$Path, {"C:\\Users\\matte\\Desktop\\Mathematica_package\\"}];
(* Import the code *)
<<"DMFT.wl";


(* Import the input file *)
<<"InputFile.wl";
(* Prepare everything *)
<<"Preparation.wl";

(*
BathParameters[[2,1]] = ConstantArray[1.0*IdentityMatrix[f], Nbath];
BathParameters[[2,2]] = ConstantArray[1.0*IdentityMatrix[f], Nbath];*)


Do[
	(* First print *)
	Print[Style["DMFT Loop n. ", 20, Bold, Red], Style[DMFTiterator, 20, Bold, Red]];
	Print[Style["\t\t Exact Diagonalization start", 16, Bold, Orange]];
	Print["Bath parameters: ", BathParameters];

	(* Build and diagonalize the AIM Hamiltonian + print timing *)
	Print["E.D. time: ", First @ AbsoluteTiming[
	
		Hsectors = HImp[HnonlocBlocks, HlocBlocks, BathParameters];
		(*checkhermiticity = Apply[And, HermitianMatrixQ[#] &/@ Hsectors];
		If[!checkhermiticity, Abort[]; ];*)
		
		eigs = Map[
			Eigs[
				Hsectors[[#]], 
				"Temperature" -> T, 
				"MinLanczosDim" -> MinLanczosDim, 
				"DegeneracyThreshold" -> DegeneracyThreshold,
				(*"MaxIterations" -> MaxLanczosIter, *)
				"MinEigenvalues" -> MinNumberOfEigs
			]&, (QnsSectorListToDiagonalize /. SectorsDispatch)
		];
		EgsSectorList = eigs[[All, 1]];
		GsSectorList = eigs[[All, 2]];
		
	]," sec.\n"];
	
	
	(* -------------------------- *)
	(* ONLY ZERO TEMPERATURE CALCULATIONS *)
	(* -------------------------- *)
	If[T != 0, Print["Error. T != 0 not supported."]; Break[]; ];
		
	(* COMPUTE THE GROUND STATE *)
	(* ground state energy (lowest of all the sectors) *)
	Egs = Min[Flatten[EgsSectorList]];
	(* sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy *)
	GsSectorIndex = Position[
		EgsSectorList,
		_?((Abs[# - Egs] < DegeneracyThreshold)&)
	];
	(* list of all the degenerate ground states *)
	Gs = Apply[GsSectorList[[##]]&, GsSectorIndex, {1}];
	(* list of quantum numbers of the degenerate ground states *)
	GsQns = QnsSectorListToDiagonalize[[GsSectorIndex[[All, 1]]]];	
		
	(* print relevant information about the ground state and observables *)
	Print["\t\t Ground state info:\n", "Egs = ", Egs, "   Quantum numbers = ", GsQns];
	Print["\n\t\t Observables:"];
	cdgc = Table[ Which[
		\[Sigma] != \[Rho],
		If[EdMode == "Magnetic", 
			0,
		(* else EdMode Raman *) 
			CdgC[L, f, Norb, {1,1}, {\[Sigma], \[Rho]}, {orb, orb}, Sectors, EgsSectorList, GsSectorList, T] 
		],
		\[Sigma] == \[Rho],
		Density[L, f, Norb, 1, \[Sigma], orb, Sectors[[QnsSectorListToDiagonalize/.SectorsDispatch]], EgsSectorList, GsSectorList, T]
	], {orb, Norb}, {\[Sigma], f}, {\[Rho], f}];
	Table[
		Print["\!\(\*SuperscriptBox[\(c\), \(\[Dagger]\)]\)_",orb,"c_",orb," = ", MatrixForm[cdgc[[orb]]]]
	, {orb, Norb}];
	density = Tr[#] &/@ cdgc;
	Print["Impurity density = ", density, "; total density = ", Total[density] ];
	docc = Table[ 
		Sum[ SquareDensity[L, f, Norb, {1,1}, {\[Sigma],\[Rho]}, {orb,orb}, Sectors[[QnsSectorListToDiagonalize/.SectorsDispatch]], EgsSectorList, GsSectorList, T], {\[Sigma], f}, {\[Rho], \[Sigma]+1, f}]
	, {orb, Norb}];
	Print["Impurity double occupancy = ", docc ];
		
	(* identify independent parameters, i.e. the minimal set of bath parameters that you need to compute stuff *)
	IndependentParameters = TakeIndependentParameters[BathParameters, independentsymbolsindexes];			
		
	(* COMPUTE GREEN FUNCTIONS *)
	Which[
	(* ----------------------------------------------------------------------------------------- *)
	(* if there is orbital symmetry, you compute many body functions just for ONE representative orbital *)
	(* ----------------------------------------------------------------------------------------- *)
		OrbitalSymmetry,
		Print["\n Green functions calculation time: ", First @ AbsoluteTiming[	
			(* G(i\[Omega]) *)
			Gimp = Mean[Apply[
				GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, MinLanczosMomenta, EdMode, i\[Omega]]&,
				{Gs, GsQns}\[Transpose]
			, {1}]];
			(* G^-1(i\[Omega]) *)
			InverseG = InverseGreenFunction[Gimp, EdMode];
			(* G_0^-1(i\[Omega]) *)
			If[DMFTiterator == 1, 
				Weiss = WeissFieldNew[L, f, Norb, H0, symbols, z, EdMode][[1]]; 
			];
			(* -- DELETE -- Gimp = 1./(Weiss/.Thread[independentsymbols -> IndependentParameters])/.{z -> #} &/@ i\[Omega];*)
			(* -- DELETE -- InverseG = InverseGreenFunction[Gimp, EdMode];*)
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0]; (* memorize the previous one *)
			InverseG0 = (Weiss/.Thread[independentsymbols -> IndependentParameters])/.{z -> #} &/@ i\[Omega];
			(* \[CapitalSigma](i\[Omega]) *)
			\[CapitalSigma]old = If[DMFTiterator == 1, 0*InverseG, \[CapitalSigma]];
			\[CapitalSigma] = InverseG0 - InverseG;
			(* G_loc(i\[Omega]) *)
			LocalGold = If[DMFTiterator == 1, 0*InverseG, LocalG];
			LocalG = LocalGreenFunction[LatticeEnergies[[All,1,1]], LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode, SublatticesQ];
			(* Weiss field (numerical) *)
			WeissNumeric = WeissFieldNumeric[
				W[[1]], \[Mu] - \[Delta][[1]], LocalG, LocalGold, \[CapitalSigma], \[CapitalSigma]old, i\[Omega], EdMode, 
				Mix -> If[DMFTiterator>2, Mixing, 0.0],
				Lattice -> LatticeType, 
				LatticeDimension -> LatticeDim
			];		
			(* compute and print z *)
			Z = QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode, FitCutoff -> 50];
			Print["Quasiparticle weight z = ", Z];			
		], " sec."];
			
			
		(* SELF CONSISTENCY *)
		Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
		Print["S.C. time: ", First @ AbsoluteTiming[	
			
		BathParameters = symbols /. SelfConsistency[
				Weiss, independentsymbols, z, IndependentParameters, WeissNumeric, i\[Omega], EdMode,
				Minimum -> MinimizationType, 
				Method -> MinimizationMethod,
				NumberOfFrequencies -> CGNMatsubara, 
				MaxIterations -> CGMaxIterations, 
				AccuracyGoal -> CGAccuracy,
				FitWeight -> CGWeight
			]
		
		], " sec." ];
			
		(* compute error *)
		error = DMFTError[InverseG0, InverseG0old, EdMode];,
		
	(* -------------------- *)	
	(* if NO ORBITAL SYMMETRY *)
	(* -------------------- *)
		!OrbitalSymmetry,
		Print["\n Green functions calculation time: ", First @ AbsoluteTiming[
		(* { G(i\[Omega])_orb=1 , G(i\[Omega])_orb=2 , ...} *)
		Gimp = Table[
			Mean[Apply[
				GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, MinLanczosMomenta, EdMode, i\[Omega]]&,
				{Gs, GsQns}\[Transpose]
			, {1}]], {orb, Norb}];
		(* G^-1(i\[Omega])_orb=1 , G^-1(i\[Omega])_orb=2 ... *)
		InverseG = InverseGreenFunction[#, EdMode] &/@ Gimp;
		(* { Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=1] , Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=2] , ...} *)
		If[DMFTiterator == 1, 
			Weiss = WeissFieldNew[L, f, Norb, H0, symbols, z, EdMode]
		];
		InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0];
		InverseG0 = Table[(
			Weiss[[orb]]/.Thread[independentsymbols -> IndependentParameters])/.{z -> #}&/@i\[Omega]
		, {orb, Norb}];
		(* { \[CapitalSigma]Subscript[(i\[Omega]), orb=1] , \[CapitalSigma]Subscript[(i\[Omega]), orb=2] , ...} *)
		\[CapitalSigma]old = If[DMFTiterator == 1, 0*InverseG, \[CapitalSigma]];
		\[CapitalSigma] = InverseG0 - InverseG;
		(* { Subscript[G, loc]Subscript[(i\[Omega]), orb=1] , Subscript[G, loc]Subscript[(i\[Omega]), orb=2] , ...} *)
		LocalGold = If[DMFTiterator == 1, 0*InverseG, LocalG];
		LocalG = Table[
			LocalGreenFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma][[orb]], i\[Omega], EdMode, SublatticesQ]
		, {orb, Norb}];
		(* Weiss field (numerical) *)
		WeissNumeric = Table[
			WeissFieldNumeric[
				W[[orb]], \[Mu] - \[Delta][[orb]], LocalG[[orb]], LocalGold[[orb]], \[CapitalSigma][[orb]], \[CapitalSigma]old[[orb]], i\[Omega], EdMode, 
				Mix -> If[DMFTiterator>2, Mixing, 0.0],
				Lattice -> LatticeType, 
				LatticeDimension -> LatticeDim
			]
		, {orb, Norb}];
			
		(* z *)
		Z = Table[ QuasiparticleWeight[\[CapitalSigma][[orb]], i\[Omega], EdMode, FitCutoff -> 50], {orb, Norb} ];
		Print["Quasiparticle weight z = ", Z];
		
		], " sec."];
			
			
		(* SELF CONSISTENCY *)
		Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
		Print["S.C. time: ", First @ AbsoluteTiming[
			
		BathParameters = symbols /. Flatten[Table[
			SelfConsistency[
				Weiss[[orb]], 
				IndependentSymbols[symbols[[All,orb,All]]], 
				z, 
				TakeIndependentParameters[BathParameters, IndependentSymbolsIndexes[symbols, IndependentSymbols[symbols[[All,orb,All]]]]], 
				WeissNumeric[[orb]], 
				i\[Omega], 
				EdMode,
				Minimum -> MinimizationType, 
				Method -> MinimizationMethod,
				NumberOfFrequencies -> CGNMatsubara, 
				MaxIterations -> CGMaxIterations, 
				AccuracyGoal -> CGAccuracy,
				FitWeight -> CGWeight
			]
		, {orb, Norb}]]
			
		], " sec." ];
		
		error = (1./Norb)*Sum[
			DMFTError[InverseG0[[orb]], InverseG0old[[orb]], EdMode],
		{orb, Norb}];
		
	];
	
	(* store new bath parameters and error *)
	WriteOutput[True, OutputDirectory, "hamiltonian_restart", BathParameters];
	AppendTo[ErrorList, error];
	WriteOutput[True, OutputDirectory, "error", ErrorList];
	Print["DMFT error: ", ScientificForm[error]];
	
	(*Print*)
	Print[Style["\t\t Self Consistency completed", 16, Bold, Magenta]];
	Print["----------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------"];

	(* Exit DMFT Loop if convengerce is reached *)
	If[DMFTiterator > DMFTMinIterations && error < DMFTerror && LastIteration, Converged = True; Break[];];
	If[error < DMFTerror, LastIteration = True, (*else*) LastIteration = False];

, {DMFTiterator, DMFTMaxIterations}]


(* POST PROCESSING *)
Which[
	OrbitalSymmetry,
	(* compute fit *)
	fitcheck = {
		SmartInverse[#] &/@ InverseG0,
		SmartInverse[#] &/@ WeissNumeric
	};
	(* <\[Psi]dg_k \[Psi]_k> *)
	\[Psi]dg\[Psi] = PsidgPsi[
		Range[LatticePoints], 
		LatticeEnergies[[All,1,1]], 
		\[Mu], \[CapitalSigma], i\[Omega]
	];
	(* flavor current *)
	flavorcurrent = Table[
		FlavorCurrent[W[[1]], \[Gamma][[1]], \[Sigma], a, \[Psi]dg\[Psi], LatticeType, LatticeDim, LatticePoints, SublatticesQ]
	, {\[Sigma], f}, {a, LatticeDim}];
	Do[
		Print["I_\[Sigma]="<>If[f==2, Which[\[Sigma]==1,"\[DownArrow]", \[Sigma]==2,"\[UpArrow]"], ToString[\[Sigma]]]<>",a="<>ToString[a]<>" = ", flavorcurrent[[\[Sigma],a]]];
	, {\[Sigma], f}, {a, LatticeDim}];
	Print[
		"(1/N)\!\(\*SubscriptBox[\(\[Sum]\), \(k\)]\) <\!\(\*SuperscriptBox[SubscriptBox[\(\[Psi]\), \(k\)], \(\[Dagger]\)]\) \!\(\*SubscriptBox[\(\[Psi]\), \(k\)]\)> = ", 
		MatrixForm[Table[Re[Mean[\[Psi]dg\[Psi][[All,\[Sigma],\[Rho]]]]], {\[Sigma], f}, {\[Rho], f}]]
	];
	(* Gimp(\[Omega]) *)
	Gimprealfreq = Mean[Apply[
		GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, MinLanczosMomenta, EdMode, \[Omega] + I*\[Eta]]&,
		{Gs, GsQns}\[Transpose]
	, {1}]];
	(* G^-1(\[Omega]) *)
	InverseGrealfreq = InverseGreenFunction[Gimprealfreq, EdMode];
	(* compute G0(\[Omega])^-1 *)
	InverseG0realfreq = (Weiss/.Thread[independentsymbols -> IndependentParameters])/.{z -> #}&/@(\[Omega] + I*\[Eta]);
	(* \[CapitalSigma](\[Omega]) *)
	\[CapitalSigma]realfreq = InverseG0realfreq - InverseGrealfreq;
	(* compute the lattice spectral function A(\[Omega]) *)
	spectralfunction = SpectralFunction[LatticeEnergies[[All, 1, 1]], LatticeWeights, \[Mu], \[CapitalSigma]realfreq, \[Omega]+I*\[Eta], EdMode, SublatticesQ];,
(* ------------------------------------------------------------------------------------ *)
	!OrbitalSymmetry,
	(* fit *)
	fitcheck = {
		Table[ SmartInverse[#] &/@ InverseG0[[orb]], {orb, Norb}],
		Table[ SmartInverse[#] &/@ WeissNumeric[[orb]], {orb, Norb}]
	};
	(* <\[Psi]dg_k \[Psi]_k> *)
(*	\[Psi]dg\[Psi] = ParallelTable[ PsidgPsi[
		Range[LatticePoints], 
		LatticeEnergies[[All,orb,orb]], 
		\[Mu], \[CapitalSigma][[orb]], i\[Omega]
	], {orb, Norb}];
	Do[
		Print[
			"(1/N)\!\(\*SubscriptBox[\(\[Sum]\), \(k\)]\) <\!\(\*SuperscriptBox[SubscriptBox[\(\[Psi]\), \(k\)], \(\[Dagger]\)]\) \!\(\*SubscriptBox[\(\[Psi]\), \(k\)]\)> (orb " <> ToString[orb] <> ") = ", 
			MatrixForm[Table[Re[Mean[\[Psi]dg\[Psi][[orb,All,\[Sigma],\[Rho]]]]], {\[Sigma], f}, {\[Rho], f}]]
		];
	, {orb, Norb}];
	(* flavor current *)
	flavorcurrent = ParallelTable[
		FlavorCurrent[W[[orb]], \[Gamma][[orb]], \[Sigma], a, \[Psi]dg\[Psi][[orb]], LatticeType, LatticeDim, LatticePoints, SublatticesQ]
	, {orb, Norb}, {\[Sigma], f}, {a, LatticeDim}];
	Do[
		Print["I_orb="<>ToString[orb]<>",\[Sigma]="<>If[f==2, Which[\[Sigma]==1,"\[DownArrow]", \[Sigma]==2,"\[UpArrow]"], ToString[\[Sigma]]]<>",a="<>ToString[a]<>" = ", flavorcurrent[[orb,\[Sigma],a]]];
	, {orb, Norb}, {\[Sigma], f}, {a, LatticeDim}];*)
	(* Gimp(\[Omega]) *)
	Gimprealfreq = Table[
		Mean[Apply[
			GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, 200, EdMode, \[Omega] + I*\[Eta]]&,
			{Gs, GsQns}\[Transpose]
		, {1}]], {orb, Norb}];
	(* compute G(\[Omega])^-1 *)
	InverseGrealfreq = InverseGreenFunction[#, EdMode] &/@ Gimprealfreq;
	(* compute G0(\[Omega])^-1 *)
	InverseG0realfreq = ParallelTable[(
		Weiss[[orb]]/.Thread[independentsymbols -> IndependentParameters])/.{z -> #}&/@(\[Omega] + I*\[Eta])
	, {orb, Norb}];
	(* compute \[CapitalSigma](\[Omega]) *)
	\[CapitalSigma]realfreq = InverseG0realfreq - InverseGrealfreq;
	(* compute the lattice spectral function A(\[Omega]) *)
	spectralfunction = ParallelTable[
		SpectralFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma]realfreq[[orb]], \[Omega]+I*\[Eta], EdMode, SublatticesQ]
	, {orb, Norb}];
]

(* print density correlations *)
Table[
	Table[
		Print[
			"<n_orb=",orb1,"_spin=",\[Sigma]," n_orb=",orb2,"_spin=",\[Rho],"> = ",
			SquareDensity[L, f, Norb, {1,1}, {\[Sigma],\[Rho]}, {orb1,orb2}, Sectors[[QnsSectorListToDiagonalize/.SectorsDispatch]], EgsSectorList, GsSectorList, T]
		]
	, {\[Sigma], f}, {\[Rho], f}]
, {orb1, Norb}, {orb2, Norb}];

(* store output *)
WriteOutput[False, OutputDirectory, "fit", fitcheck];
WriteOutput[False, OutputDirectory, "density", density];
WriteOutput[False, OutputDirectory, "double_occupancy", docc];
WriteOutput[False, OutputDirectory, "flavor_current", flavorcurrent];
WriteOutput[False, OutputDirectory, "quasiparticle_weight", Z];
WriteOutput[False, OutputDirectory, "self_energy", \[CapitalSigma]];
WriteOutput[False, OutputDirectory, "self_energy_real_frequency", \[CapitalSigma]realfreq];
WriteOutput[False, OutputDirectory, "spectral_function", spectralfunction];
WriteOutput[False, OutputDirectory, "psidgpsi", \[Psi]dg\[Psi]];


ListPlot[{
	Im[ fitcheck[[1, 2, All, 1]] ],
	Im[ fitcheck[[2, 2, All, 1]] ]
}, Joined->True, PlotMarkers->Automatic, PlotStyle->{Thick,Dashed}, PlotRange->{{0,50}, All}, AxesLabel->{"n"," "}]


PlotSpectralFunction[spectralfunction[[1]]]



(*
Table[
	CdgC[L, f, Norb, {1,1}, {\[Sigma],\[Rho]}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T]
, {orb, Norb}, {\[Sigma],f}, {\[Rho],f}]

< n1 n2 > = 
< (n1up + n1dw)(n2up + n2dw) >=
< n1up n2up > + < n1up n2dw > + < n1dw n2up > + < n1dw n2dw >
*)



ListPlot[{
Partition[
	\[Psi]dg\[Psi][[1,All,1,1]] 
	, Sqrt[LatticePoints]
]\[Transpose][[1,All]],
Partition[
	\[Psi]dg\[Psi][[2,All,1,1]] 
	, Sqrt[LatticePoints]
]\[Transpose][[1,All]]
},
DataRange->{-Pi,Pi},
PlotRange->{0,1}, Joined->True, 
Filling->Axis,
PlotLabel->"s=5, \[CapitalOmega] = 0.5D, \[Gamma]=\[Pi]/2",
Axes->False, Frame->True,
FrameStyle->Directive[Black,16],
AspectRatio->0.5,
Epilog->Line[{{0,0},{0,1}}],
FrameLabel->{"\!\(\*SubscriptBox[\(k\), \(x\)]\)a",None},
PlotLegends->{"\!\(\*SubscriptBox[\(n\), \(g\)]\)(\!\(\*SubscriptBox[\(k\), \(x\)]\),\!\(\*SubscriptBox[\(k\), \(y\)]\)=0)","\!\(\*SubscriptBox[\(n\), \(e\)]\)(\!\(\*SubscriptBox[\(k\), \(x\)]\),\!\(\*SubscriptBox[\(k\), \(y\)]\)=0)"}

]


QnsSectorListToDiagonalize = Select[
	QnsSectorList,
	(
		#[[1]] + #[[2]] == 5
		&& ((#[[3]] == 2 && #[[4]] == 3) || (#[[3]] == 3 && #[[4]] == 2))
	)&
];

QnsSectorList[[
	QnsSectorListToDiagonalize /. SectorsDispatch
]] == QnsSectorListToDiagonalize
