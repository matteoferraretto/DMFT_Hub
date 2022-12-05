(* ::Package:: *)

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
(* 1. compute G^-1 in real frequencies *)
Gimprealfreq = Table[
	Mean[MapApply[
		GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega] + I*\[Eta]]&,
		{Gs, GsQns}\[Transpose]
	]], {orb, Norb}];
InverseGrealfreq = InverseGreenFunction[Gimprealfreq, EdMode];
(* 2. compute Subscript[G, 0]^-1 in real frequencies *)
InverseG0realfreq = Table[(
	(Weiss/.{\[Mu]eff -> \[Mu] - \[Delta][[orb]]})/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode]])/.{z -> #}&/@(\[Omega] + I*\[Eta])
, {orb, Norb}];
(* 3. compute \[CapitalSigma] in real frequencies *)
\[CapitalSigma]realfreq = InverseG0realfreq - InverseGrealfreq;
(* 4. compute the lattice spectral function *)
spectralfunction = Table[
	SpectralFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma]realfreq[[orb]], \[Omega]+I*\[Eta], EdMode]
, {orb, Norb}];


(* plot stuff *)
Which[
	LatticeType == "Bethe",
	(* spectral function *)
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
	(* energy-resolved spectral function A(\[Epsilon], \[Omega]) *)
	(* *)
	spectralfunctionresolved = MomentumResolvedSpectralFunction[
		LatticeEnergies[[All, 1, 1]][[1;;-9;;10]], \[Mu], \[CapitalSigma]realfreq[[1]][[1;;-1;;10]], Range[Length[LatticeEnergies[[All, 1, 1]][[1;;-9;;10]]]], (\[Omega]+I*\[Eta])[[1;;-1;;10]], EdMode
	];
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
	]
]
(* check if the integral gives 1 *)
Print["Integral of spectral function = ", 
	d\[Omega] * Total[#] &/@ Table[spectralfunction[[orb]][[All,2]] , {orb, Norb}]
];


Range[LatticePoints][[;;-9;;10]]
