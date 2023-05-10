(* ::Package:: *)

(* Compute many body functions for all the Matsubara frequencies *)
(* get all Matsubara frequencies *)
i\[Omega] = Table[(2n-1)Pi*I*TMats, {n, NMatsubara}];
(* initialize real frequencies *)
\[Omega] = Table[\[Omega]min + n*d\[Omega], {n, 0, NReal}];

(* compute converged many body functions in Matsubara and real frequencies *)
Which[
	(EdMode == "Normal" || EdMode == "Superc" || EdMode == "Raman") && OrbitalSymmetry && !SublatticesQ,
	(* G(i\[Omega]) *)
	Gimp = Mean[Apply[
		GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, MinLanczosMomenta, EdMode, i\[Omega]]&,
		{Gs, GsQns}\[Transpose]
	, {1}]];
	(* compute G(\[Omega]) for all orbitals *)
	Gimprealfreq = Mean[Apply[
		GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, MinLanczosMomenta, EdMode, \[Omega] + I*\[Eta]]&,
		{Gs, GsQns}\[Transpose]
	, {1}]];
	(* G^-1(i\[Omega]) *)
	InverseG = InverseGreenFunction[Gimp, EdMode];
	(* G^-1(\[Omega]) *)
	InverseGrealfreq = InverseGreenFunction[Gimprealfreq, EdMode];
	(* G_0^-1(i\[Omega]) *)
	InverseG0 = (Weiss/.Thread[symbols -> IndependentParameters])/.{z -> #}&/@i\[Omega];
	(* compute G0(\[Omega])^-1 *)
	InverseG0realfreq = (Weiss/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode]])/.{z -> #}&/@(\[Omega] + I*\[Eta]);
	(* \[CapitalSigma](i\[Omega]) *)
	\[CapitalSigma] = InverseG0 - InverseG;
	(* \[CapitalSigma](\[Omega]) *)
	\[CapitalSigma]realfreq = InverseG0realfreq - InverseGrealfreq;
	(* G_loc(i\[Omega]) *)
	LocalG = LocalGreenFunction[LatticeEnergies[[All,1,1]], LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode, SublatticesQ];
	(* compute the lattice spectral function A(\[Omega]) *)
	(*spectralfunction = SpectralFunction[LatticeEnergies[[All, 1, 1]], LatticeWeights, \[Mu], \[CapitalSigma]realfreq, \[Omega]+I*\[Eta], EdMode];*),
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
	(EdMode == "Normal" || EdMode == "Superc" || EdMode == "Raman") && !OrbitalSymmetry && !SublatticesQ,
	(* compute G(i\[Omega]) for all orbitals *)
	Gimp = Table[
		Mean[Apply[
			GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, MinLanczosMomenta, EdMode, i\[Omega]]&,
			{Gs, GsQns}\[Transpose]
		, {1}]], {orb, Norb}];
	(* compute G(\[Omega]) for all orbitals *)
	Gimprealfreq = Table[
		Mean[Apply[
			GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, 200, EdMode, \[Omega] + I*\[Eta]]&,
			{Gs, GsQns}\[Transpose]
		, {1}]], {orb, Norb}];
	(* {G^-1(i\[Omega])_orb=1 , G^-1(i\[Omega])_orb=2 ...} *)
	InverseG = InverseGreenFunction[#, EdMode] &/@ Gimp;
	(* compute G(\[Omega])^-1 *)
	InverseGrealfreq = InverseGreenFunction[#, EdMode] &/@ Gimprealfreq;
	(* {G0^-1(i\[Omega])_orb=1, G0^-1(i\[Omega])_orb=2 ...} *)
	InverseG0 = Table[(
		Weiss[[orb]]/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode]])/.{z -> #}&/@i\[Omega]
	, {orb, Norb}];
	(* compute G0(\[Omega])^-1 *)
	InverseG0realfreq = Table[(
		Weiss[[orb]]/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode]])/.{z -> #}&/@(\[Omega] + I*\[Eta])
	, {orb, Norb}];
	(* { \[CapitalSigma](i\[Omega])_orb=1 , \[CapitalSigma](i\[Omega])_orb=2 , ...} *)
	\[CapitalSigma] = InverseG0 - InverseG;
	(* compute \[CapitalSigma](\[Omega]) *)
	\[CapitalSigma]realfreq = InverseG0realfreq - InverseGrealfreq;
	(* { Gloc(i\[Omega])_orb=1 , Gloc(i\[Omega])_orb=2 , ...} *)
	LocalG = Table[
		LocalGreenFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma][[orb]], i\[Omega], EdMode, SublatticesQ]
	, {orb, Norb}];
	(* compute the lattice spectral function A(\[Omega]) *)
	spectralfunction = Table[
		SpectralFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma]realfreq[[orb]], \[Omega]+I*\[Eta], EdMode]
	, {orb, Norb}];,
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
	(EdMode == "InterorbSuperc" || EdMode == "FullSuperc") && !OrbitalSymmetry && !SublatticesQ,
	Gimp = Mean[Apply[
		GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, MinLanczosMomenta, EdMode, i\[Omega]]&,
		{Gs, GsQns}\[Transpose]
	, {1}]];
	InverseG = InverseGreenFunction[Gimp, EdMode];
	(* G0^-1(i\[Omega]) *)
	InverseG0 = (
		Weiss/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode]]
	)/.{z -> #} &/@ i\[Omega];
	(* \[CapitalSigma](i\[Omega]) *)
	\[CapitalSigma] = InverseG0 - InverseG;
	(* Gloc(i\[Omega]) *)
	LocalG = LocalGreenFunction[LatticeEnergies, LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode, SublatticesQ];,
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
	(EdMode == "Normal" || EdMode == "Superc" || EdMode == "Raman") && SublatticesQ,
	(* repeat one more ED step *)
	Do[
		Hsectors = HImp[Norb, HnonlocBlocks, HlocBlocks, BathParameters[[sublattice]], InteractionParameters[[sublattice]], EdMode];
			eigs = Map[
				Eigs[#, 
					"Temperature" -> T, 
					"MinLanczosDim" -> MinLanczosDim, 
					"DegeneracyThreshold" -> DegeneracyThreshold,
					"MaxIterations" -> MaxLanczosIter, 
					"MinEigenvalues" -> MinNumberOfEigs
				]&, Hsectors];
		EgsSectorList = eigs[[All, 1]];
		GsSectorList = eigs[[All, 2]];
		ClearAll[eigs];
		Egs = Min[Flatten[EgsSectorList]];
		GsSectorIndex = Position[EgsSectorList, _?((Abs[# - Egs] < DegeneracyThreshold)&)];
		Gs = Apply[GsSectorList[[##]]&, GsSectorIndex, {1}];
		GsQns = QnsSectorList[[GsSectorIndex[[All, 1]]]];	
		(* double occupancy *)
		docc = Table[ SquareDensity[L, f, Norb, {1,1}, {1,2}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T], {orb, Norb}]
		WriteOutput[True, OutputDirectory, "double_occupancy_sublattice_"<>If[sublattice==1,"A","B"], docc];
		(* < cdg_\[Sigma] c_\[Rho] > *)
		If[EdMode == "Raman",
			cdgc = Table[
				Which[
					\[Sigma] != \[Rho], CdgC[L, f, Norb, {1,1}, {\[Sigma], \[Rho]}, {orb, orb}, Sectors, EgsSectorList, GsSectorList, T],
					\[Sigma] == \[Rho], Density[L, f, Norb, 1, \[Sigma], orb, Sectors, EgsSectorList, GsSectorList, T]
				]
			, {orb, Norb}, {\[Sigma], f}, {\[Rho], f}];
			WriteOutput[True, OutputDirectory, "cdgc_sublattice_"<>If[sublattice==1,"A","B"], cdgc];
		];
		If[OrbitalSymmetry,
			IndependentParameters[[sublattice]] = TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters[[sublattice]], EdMode];
			Gimp = Mean[Apply[
				GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, MinLanczosMomenta, EdMode, i\[Omega]]&,
				{Gs, GsQns}\[Transpose]
			, {1}]];
			InverseG = InverseGreenFunction[Gimp, EdMode];
			InverseG0 = (Weiss[[sublattice]]/.Thread[symbols -> IndependentParameters[[sublattice]]])/.{z -> #} &/@ i\[Omega];
			\[CapitalSigma][[sublattice]] = InverseG0 - InverseG;,
		(* else, not OrbitalSymmetry *)
			Do[ IndependentParameters[[orb, sublattice]] = TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters[[sublattice]], EdMode], {orb, Norb}];
			Gimp = Table[
				Mean[Apply[
					GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, MinLanczosMomenta, EdMode, i\[Omega]]&,
					{Gs, GsQns}\[Transpose]
				, {1}]], {orb, Norb}];
			InverseG = InverseGreenFunction[#, EdMode] &/@ Gimp;
			InverseG0 = Table[(
				Weiss[[orb, sublattice]]/.Thread[symbols -> IndependentParameters[[orb, sublattice]]])/.{z -> #}&/@i\[Omega]
			, {orb, Norb}];
			Do[\[CapitalSigma][[orb, sublattice]] = InverseG0[[orb]] - InverseG[[orb]];, {orb, Norb}];
		];
		(* real frequency many body functions *)
		(*Gimprealfreq = Table[
			Mean[Apply[
				GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, 200, EdMode, \[Omega] + I*\[Eta]]&,
				{Gs, GsQns}\[Transpose]
			, {1}]]
		, {orb, Norb}];
		InverseGrealfreq = InverseGreenFunction[#, EdMode] &/@ Gimprealfreq;
		InverseG0realfreq = Table[(
			Weiss[[orb, sublattice]]/.Thread[symbols -> IndependentParameters[[orb, sublattice]]])/.{z -> #}&/@(\[Omega]+I*\[Eta])
		, {orb, Norb}];
		\[CapitalSigma]realfreq = InverseG0realfreq - InverseGrealfreq;
		WriteOutput[True, OutputDirectory, "self_energy_real_frequency_sublattice_"<>If[sublattice==1,"A","B"], \[CapitalSigma]realfreq];*)
	, {sublattice, 1, 2}];
	(* Save the self energy as a Norb x 2 x NMatsubara x f x f tensor *)
	Do[ WriteOutput[True, OutputDirectory, "self_energy_sublattice_"<>If[sublattice==1,"A","B"], \[CapitalSigma][[sublattice]]], {sublattice, 1, 2}];
	(* reshape self energy to make it a Norb x NMatsubara x 2f x 2f tensor *)
	If[OrbitalSymmetry,
		\[CapitalSigma]matrix = (ArrayFlatten[{{#,0},{0,0*#}}] &/@ \[CapitalSigma][[1]]) + (ArrayFlatten[{{0*#,0},{0,#}}] &/@ \[CapitalSigma][[2]]);,
	(* else Norb > 1 *)
		\[CapitalSigma]matrix = Table[
			(ArrayFlatten[{{#,0},{0,0*#}}] &/@ \[CapitalSigma][[orb, 1]]) + (ArrayFlatten[{{0*#,0},{0,#}}] &/@ \[CapitalSigma][[orb, 2]]);
		, {orb, Norb}]
	];
	ClearAll[\[CapitalSigma]]; \[CapitalSigma] = \[CapitalSigma]matrix; ClearAll[\[CapitalSigma]matrix];
];


(* store converged observables and many body functions *)
If[!SublatticesQ,
	WriteOutput[True, OutputDirectory, "density", density];
	WriteOutput[True, OutputDirectory, "cdgc", cdgc];
	WriteOutput[True, OutputDirectory, "double_occupancy", docc];
	WriteOutput[True, OutputDirectory, "quasiparticle_weight", Z];
	If[EdMode == "Superc" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		WriteOutput[True, OutputDirectory, "intraorbital_order_parameter", \[Phi]];
		If[EdMode != "Superc", 
			WriteOutput[True, OutputDirectory, "interorbital_order_parameter", \[CapitalXi]]
		];
	];
	WriteOutput[True, OutputDirectory, "self_energy", \[CapitalSigma]];
	WriteOutput[True, OutputDirectory, "self_energy_real_frequency", \[CapitalSigma]realfreq];
	WriteOutput[True, OutputDirectory, "spectral_function", spectralfunction];
];


(* momentum resolved spectral function *)
(*If[
	OrbitalSymmetry || Norb == 1,
	spectralfunctionresolved = MomentumResolvedSpectralFunction[
		LatticeEnergies[[All, 1, 1]], 
		\[Mu] - \[Delta][[1]], 
		\[CapitalSigma]realfreq,
		HighSymmetryPath[ Length[LatticeEnergies[[All,1,1]]], LatticeType, LatticeDim],
		\[Omega] + I*\[Eta], 
		EdMode
	],
(* else, if no orbital symmetry *)
	spectralfunctionresolved = Table[
		MomentumResolvedSpectralFunction[
			LatticeEnergies[[All, orb, orb]], 
			\[Mu] - \[Delta][[orb]], 
			\[CapitalSigma]realfreq[[orb]],
			HighSymmetryPath[ Length[LatticeEnergies[[All, orb, orb]]], LatticeType, LatticeDim],
			\[Omega] + I*\[Eta], 
			EdMode
		]
	, {orb, Norb}]
];

WriteOutput[True, OutputDirectory, "momentum_resolved_spectral_function", spectralfunctionresolved];
ClearAll[spectralfunctionresolved];*)


(* specific observables of Raman coupled systems *)
(* this code should work regardless of the value of SublatticeQ, by a miracle *)
If[ EdMode == "Raman",
	If[OrbitalSymmetry,
		flavordistribution = MomentumDistributedDensityRaman[
			HighSymmetryPath[ Length[ LatticeEnergies[[All, 1, 1]] ], LatticeType, LatticeDim], 
			LatticeEnergies[[All, 1, 1]], \[Mu], \[CapitalSigma], i\[Omega]
		];
		Print["Dimensions of the tensor <cdg_k\[Alpha] c_k\[Beta]>: ", Dimensions[flavordistribution]];
		(* flavor current *)
		flavorcurrent = Table[
			FlavorCurrent[W[[1]], \[Gamma][[1]], \[Sigma], a, flavordistribution, LatticeType, LatticeDim, LatticePoints, SublatticesQ]
		, {\[Sigma], f}, {a, LatticeDim}];
		Do[
			Print["I_\[Sigma]="<>If[f==2, Which[\[Sigma]==1,"\[DownArrow]", \[Sigma]==2,"\[UpArrow]"], ToString[\[Sigma]]]<>",a="<>ToString[a]<>" = ", flavorcurrent[[\[Sigma],a]]];
		, {\[Sigma], f}, {a, LatticeDim}];
		(* kinetic energy *)
		Ekin = KineticEnergy[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega], EdMode, "OrbitalSymmetry" -> OrbitalSymmetry, "FlavorDistribution" -> {flavordistribution}];
		Print["Ekin = ", Ekin];,
	(* -------- else, if no orbital symmetry ------------- *)
		flavordistribution = Table[
			MomentumDistributedDensityRaman[
				HighSymmetryPath[ Length[ LatticeEnergies[[All, orb, orb]] ], LatticeType, LatticeDim], 
				LatticeEnergies[[All, orb, orb]], \[Mu], \[CapitalSigma][[orb]], i\[Omega]
			]
		, {orb, Norb}];
		(* flavor current *)
		flavorcurrent = Table[
			FlavorCurrent[W[[orb]], \[Gamma][[orb]], \[Sigma], a, flavordistribution[[orb]], LatticeType, LatticeDim, LatticePoints, SublatticesQ]
		, {orb, Norb}, {\[Sigma], f}, {a, LatticeDim}];
		Do[
			Print["I_orb="<>ToString[orb]<>",\[Sigma]="<>If[f==2, Which[\[Sigma]==1,"\[DownArrow]", \[Sigma]==2,"\[UpArrow]"], ToString[\[Sigma]]]<>",a="<>ToString[a]<>" = ", flavorcurrent[[orb,\[Sigma],a]]];
		, {orb, Norb}, {\[Sigma], f}, {a, LatticeDim}];
		(* kinetic energy *)
		Ekin = KineticEnergy[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega], EdMode, "FlavorDistribution" -> flavordistribution];
		Print["Ekin = ", Ekin];
	];

	WriteOutput[True, OutputDirectory, "flavor_current", flavorcurrent];
	(*WriteOutput[True, OutputDirectory, "kinetic_energy", Ekin];*)
	WriteOutput[True, OutputDirectory, "flavor_resolved_momentum_distributed_cdgc", flavordistribution];
	(*ClearAll[flavordistribution];*)
];


(* SPECTRAL FUNCTIONS *)
(*
Which[
	LatticeType == "Bethe",
	(* spectral function A(\[Omega]) *)
	Print @ Show[
		ListPlot[
			spectralfunction,
			Joined -> True, 
			PlotRange -> All,
			Filling -> Axis,
			AxesLabel -> {"\[Omega]", "DoS"},
			AxesStyle -> Directive[Black, 14],
			AspectRatio -> 1/2
		],
		ListPlot[
			Table[
				{LatticeEnergies[[All, orb, orb]], LatticePoints * LatticeWeights / (2. * W[[orb]])}\[Transpose]
			, {orb, Norb}], 
			Joined -> True, 
			PlotStyle -> Dashing[.05],
			PlotRange -> All,
			Filling -> Axis,
			AxesLabel -> {"\[Omega]", "DoS"},
			AxesStyle -> Directive[Black, 14],
			AspectRatio -> 1/2
		]
	];
	(* compute energy-resolved spectral function A(\[Epsilon], \[Omega]) *)
	Module[{
		energies = LatticeEnergies[[All, 1, 1]], (* take one element every 10 *)
		path,
		\[CapitalSigma] = If[OrbitalSymmetry,
			\[CapitalSigma]realfreq[[1;;-1;;20]],
		(* else, if no orbital sym. *)
			\[CapitalSigma]realfreq[[All, 1;;-1;;20]]
		], (* again one element out of 10 *)
		zlist = (\[Omega]+I*\[Eta])[[1;;-1;;20]]
		},
		path = HighSymmetryPath[Length[energies], LatticeType, LatticeDim];
		If[OrbitalSymmetry,
			spectralfunctionresolved = MomentumResolvedSpectralFunction[energies, \[Mu], \[CapitalSigma], path, zlist, EdMode];
			(* plot it *)
			Print @ ListDensityPlot[
				spectralfunctionresolved,
				FrameLabel -> {"\[Epsilon]", "\[Omega]"},
				FrameStyle -> Directive[Black, 14],
				PlotLegends -> Automatic,
				DataRange -> {{-W[[1]],W[[1]]}, {\[Omega]min, \[Omega]max}}
			],
		(* else, if no orbital sym. *)
			spectralfunctionresolved = Table[
				MomentumResolvedSpectralFunction[energies, \[Mu], \[CapitalSigma][[orb]], path, zlist, EdMode];
			, {orb, Norb}];
			Do[
				(* plot it *)
				Print @ ListDensityPlot[
					spectralfunctionresolved[[orb]],
					FrameLabel -> {"\[Epsilon]", "\[Omega]"},
					FrameStyle -> Directive[Black, 14],
					PlotLegends -> Automatic,
					DataRange -> {{-W[[1]],W[[1]]}, {\[Omega]min, \[Omega]max}}
				]
			, {orb, Norb}]
		];
	];,
(* ---------------------------------------- *)
(* ---------------------------------------- *)
(* ---------------------------------------- *)
	LatticeType == "Hypercubic",
	Print @ Show[
		ListPlot[
			spectralfunction,
			PlotRange -> All,
			Joined -> True,
			Filling -> Axis
		],
		Histogram[
			Table[LatticeEnergies[[All, orb, orb]], {orb, Norb}],
			{.15}, "PDF",
			ChartStyle -> Directive[Opacity[.3]]
		]
	];
	(* compute energy-resolved spectral function A(\[Epsilon], \[Omega]) *)
	spectralfunctionresolved = 
	MomentumResolvedSpectralFunction[
			LatticeEnergies[[All, 1, 1]], 
			\[Mu], 
			\[CapitalSigma]realfreq[[1]][[1;;-1;;100]], 
			HighSymmetryPath[Length[ LatticeEnergies[[All,1,1]] ], LatticeType, LatticeDim],
			\[Omega][[1;;-1;;100]] + I*\[Eta], 
			EdMode
	]; 
	(* plot it *)
	Print @ ListDensityPlot[
		spectralfunctionresolved,
		FrameLabel -> {"high symmetry path", "\[Omega]"},
		FrameStyle -> Directive[Black, 14],
		PlotLegends -> Automatic,
		DataRange -> {{0, 1}, {\[Omega]min, \[Omega]max}}
	]
]

(* check if the integral gives 1 *)
Print["Integral of spectral function = ", 
	If[OrbitalSymmetry,
		d\[Omega] * Total[ spectralfunction[[All,2]] ],
	(* else *)
		d\[Omega] * Total[#] &/@ Table[spectralfunction[[orb]][[All,2]] , {orb, Norb}]
	]
];
*)
