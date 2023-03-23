(* ::Package:: *)

(* Plot DMFT error *)
Print @ ListLogPlot[
	ErrorList,
	Axes->False, Frame->True,
	FrameLabel -> {"iteration", "DMFT error"},
	FrameStyle -> Directive[Black,14],
	Joined -> True, PlotMarkers -> Automatic,
	FrameTicks -> Automatic,
	PlotStyle -> PointSize[.025]
]


(* Plot momentum-resolved spectral function *)
(* spectralfunctionresolved = Import[OutputDirectory<>"momentum_resolved_spectral_function.m"];
If[
	OrbitalSymmetry || Norb == 1,
	Print @ PlotSpectralFunctionRaman[Abs[spectralfunctionresolved], \[Omega], \[Mu], LatticeType, LatticeDim],
(* else if no orbital symmetry *)
	Do[
		Print @ PlotSpectralFunctionRaman[Abs[spectralfunctionresolved[[orb]]], \[Omega], \[Mu], LatticeType, LatticeDim]
	, {orb, Norb}]	
]; *)


(* plot the spectral function *)
spectralfunction = Import[OutputDirectory<>"spectral_function.m"];
If[
	OrbitalSymmetry || Norb == 1,
	Print @ ListPlot[
		spectralfunction,
		Joined -> True, 
		PlotRange -> All,
		Filling -> Axis,
		AxesLabel -> {"\[Omega]", "DoS"},
		AxesStyle -> Directive[Black, 14],
		AspectRatio -> 1/2
	],
(* else if no orbital symmetry *)
	Do[
		Print @ ListPlot[
			spectralfunction[[orb]],
			Joined -> True, 
			PlotRange -> All,
			Filling -> Axis,
			AxesLabel -> {"\[Omega]", "DoS"},
			AxesStyle -> Directive[Black, 14],
			AspectRatio -> 1/2
		]
	, {orb, Norb}]
]
