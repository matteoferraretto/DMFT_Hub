(* ::Package:: *)

(* Compute many body functions for all the Matsubara frequencies *)
(* get all Matsubara frequencies *)
i\[Omega] = Table[(2n-1)Pi*I*TMats, {n, NMatsubara}];
(* initialize real frequencies *)
\[Omega] = Table[\[Omega]min + n*d\[Omega], {n, 0, NReal}];

(* save observables on file *)
WriteOutput[Converged, OutputDirectory, "density", density];
WriteOutput[Converged, OutputDirectory, "double_occupancy", docc];
WriteOutput[Converged, OutputDirectory, "quasiparticle_weight", Z];
If[EdMode == "Superc" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
	WriteOutput[Converged, OutputDirectory, "intraorbital_order_parameter", \[Phi]];
	If[EdMode != "Superc", 
		WriteOutput[Converged, OutputDirectory, "interorbital_order_parameter", \[CapitalXi]]
	];
];

(* compute converged many body functions *)
Which[
	(EdMode == "Normal" || EdMode == "Superc") && OrbitalSymmetry,
	(* G(i\[Omega]) *)
	Gimp = Mean[Apply[
		GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
		{Gs, GsQns}\[Transpose]
	, {1}]];
	(* compute G(\[Omega]) for all orbitals *)
	Gimprealfreq = Mean[Apply[
		GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega] + I*\[Eta]]&,
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
	LocalG = LocalGreenFunction[LatticeEnergies[[All,1,1]], LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode];
	(* compute the lattice spectral function A(\[Omega]) *)
	spectralfunction = SpectralFunction[LatticeEnergies[[All, 1, 1]], LatticeWeights, \[Mu], \[CapitalSigma]realfreq, \[Omega]+I*\[Eta], EdMode];,
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
	(EdMode == "Normal" || EdMode == "Superc") && !OrbitalSymmetry,
	(* compute G(i\[Omega]) for all orbitals *)
	Gimp = Table[
		Mean[Apply[
			GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
			{Gs, GsQns}\[Transpose]
		, {1}]], {orb, Norb}];
	(* compute G(\[Omega]) for all orbitals *)
	Gimprealfreq = Table[
		Mean[Apply[
			GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega] + I*\[Eta]]&,
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
		LocalGreenFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma][[orb]], i\[Omega], EdMode]
	, {orb, Norb}];
	(* compute the lattice spectral function A(\[Omega]) *)
	spectralfunction = Table[
		SpectralFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma]realfreq[[orb]], \[Omega]+I*\[Eta], EdMode]
	, {orb, Norb}];,
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
	(EdMode == "InterorbSuperc" || EdMode == "FullSuperc") && !OrbitalSymmetry,
	Gimp = Mean[Apply[
		GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
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
	LocalG = LocalGreenFunction[LatticeEnergies, LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode];
];

(* store converged many body functions *)
WriteOutput[True, OutputDirectory, "self_energy", \[CapitalSigma]];
WriteOutput[True, OutputDirectory, "self_energy_real_frequency", \[CapitalSigma]realfreq];


(* Plot DMFT error *)
(*
Print @ ListLogPlot[
	ErrorList,
	Axes->False, Frame->True,
	FrameLabel->{"iteration", "DMFT error"},
	FrameStyle->Directive[Black,14],
	Joined->True, PlotMarkers->Automatic,
	PlotStyle->PointSize[.025]
]
*)


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
