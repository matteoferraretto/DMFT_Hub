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



