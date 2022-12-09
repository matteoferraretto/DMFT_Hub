(* ::Package:: *)

(* Compute many body functions for all the Matsubara frequencies *)
(* get all Matsubara frequencies *)
i\[Omega] = Table[(2n-1)Pi*I*TMats, {n, NMatsubara}];

Which[
	(EdMode == "Normal" || EdMode == "Superc") || OrbitalSymmetry,
	(* G(i\[Omega]) *)
	Gimp = Mean[MapApply[
		GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
		{Gs, GsQns}\[Transpose]
	]];
	(* G^-1(i\[Omega]) *)
	InverseG = InverseGreenFunction[Gimp, EdMode];
	(* G_0^-1(i\[Omega]) *)
	InverseG0 = ((Weiss/.{\[Mu]eff -> \[Mu] - \[Delta][[1]]})/.Thread[symbols -> IndependentParameters])/.{z -> #}&/@i\[Omega];
	(* \[CapitalSigma](i\[Omega]) *)
	\[CapitalSigma] = InverseG0 - InverseG;
	(* G_loc(i\[Omega]) *)
	LocalG = LocalGreenFunction[LatticeEnergies[[All,1,1]], LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode];,
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
	(EdMode == "Normal" || EdMode == "Superc") && !OrbitalSymmetry,
	(* {G_orb=1, G_orb=2, ...} *)
	Gimp = Table[
		Mean[MapApply[
			GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
			{Gs, GsQns}\[Transpose]
		]], {orb, Norb}];
	(* {G^-1(i\[Omega])_orb=1 , G^-1(i\[Omega])_orb=2 ...} *)
	InverseG = InverseGreenFunction[#, EdMode] &/@ Gimp;
	(* {G0^-1(i\[Omega])_orb=1, G0^-1(i\[Omega])_orb=2 ...} *)
	InverseG0 = Table[(
		(Weiss/.{\[Mu]eff -> \[Mu] - \[Delta][[orb]]})/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode]])/.{z -> #}&/@i\[Omega]
	, {orb, Norb}];
	(* { \[CapitalSigma](i\[Omega])_orb=1 , \[CapitalSigma](i\[Omega])_orb=2 , ...} *)
	\[CapitalSigma] = InverseG0 - InverseG;
	(* { Gloc(i\[Omega])_orb=1 , Gloc(i\[Omega])_orb=2 , ...} *)
	LocalG = Table[
		LocalGreenFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma][[orb]], i\[Omega], EdMode]
	, {orb, Norb}];,
(* ------------------------------------------------------------------------------------ *)
(* ------------------------------------------------------------------------------------ *)
	(EdMode == "InterorbSuperc" || EdMode == "FullSuperc") && !OrbitalSymmetry,
	Gimp = Mean[MapApply[
		GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
		{Gs, GsQns}\[Transpose]
	]];
	InverseG = InverseGreenFunction[Gimp, EdMode];
	(* G0^-1(i\[Omega]) *)
	InverseG0 = (
		(Weiss/.{\[Mu]eff -> \[Mu] - \[Delta]})/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode]]
	)/.{z -> #} &/@ i\[Omega];
	(* \[CapitalSigma](i\[Omega]) *)
	\[CapitalSigma] = InverseG0 - InverseG;
	(* Gloc(i\[Omega]) *)
	LocalG = LocalGreenFunction[LatticeEnergies, LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode];
]


(* Plot DMFT error *)
Print @ ListLogPlot[
	ErrorList,
	Axes->False, Frame->True,
	FrameLabel->{"iteration", "DMFT error"},
	FrameStyle->Directive[Black,14],
	Joined->True, PlotMarkers->Automatic,
	PlotStyle->PointSize[.025]
]


(* Compute, save and plot spectral function *)
(* initialize real frequencies *)
\[Omega] = Table[\[Omega]min + n*d\[Omega], {n, 0, NReal}];

(* 1. compute G(\[Omega]) *)
Gimprealfreq = Table[
	Mean[MapApply[
		GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega] + I*\[Eta]]&,
		{Gs, GsQns}\[Transpose]
	]], {orb, Norb}];
(* compute G(\[Omega])^-1 *)
InverseGrealfreq = InverseGreenFunction[#, EdMode] &/@ Gimprealfreq;
(* 2. compute G0(\[Omega])^-1 *)
InverseG0realfreq = Table[(
	(Weiss/.{\[Mu]eff -> \[Mu] - \[Delta][[orb]]})/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode]])/.{z -> #}&/@(\[Omega] + I*\[Eta])
, {orb, Norb}];
(* 3. compute \[CapitalSigma](\[Omega]) *)
\[CapitalSigma]realfreq = InverseG0realfreq - InverseGrealfreq;
(* 4. compute the lattice spectral function *)
spectralfunction = Table[
	SpectralFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma]realfreq[[orb]], \[Omega]+I*\[Eta], EdMode]
, {orb, Norb}];


(* plot stuff *)
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
		\[CapitalSigma] = \[CapitalSigma]realfreq[[1]][[1;;-1;;20]], (* again one element out of 10 *)
		zlist = (\[Omega]+I*\[Eta])[[1;;-1;;20]]
		},
		path = HighSymmetryPath[Length[energies], LatticeType, LatticeDim];
		spectralfunctionresolved = MomentumResolvedSpectralFunction[energies, \[Mu], \[CapitalSigma], path, zlist, EdMode];
	];
	(* plot it *)
	Print @ ListDensityPlot[
		spectralfunctionresolved,
		FrameLabel -> {"\[Epsilon]", "\[Omega]"},
		FrameStyle -> Directive[Black, 14],
		PlotLegends -> Automatic,
		DataRange -> {{-W[[1]],W[[1]]}, {\[Omega]min, \[Omega]max}}
	],
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
	d\[Omega] * Total[#] &/@ Table[spectralfunction[[orb]][[All,2]] , {orb, Norb}]
];


Dimensions[Gimprealfreq]
Dimensions[InverseGrealfreq]



