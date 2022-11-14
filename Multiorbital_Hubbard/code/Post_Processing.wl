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
spectralfunction = SpectralFunction[L, f, Norb, 1, 1, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega], \[Eta]];
Which[
	LatticeType == "Bethe",
	Print @ ListPlot[
		{
			{LatticeEnergies[[All,1,1]], DoSBethe[LatticeEnergies[[All,1,1]], W[[1]]]}\[Transpose],
			spectralfunction
		}, 
		Joined->True, PlotRange->All,
		Filling->Axis,
		AxesLabel->{"\[Omega]", "DoS"},
		AxesStyle->Directive[Black, 14],
		AspectRatio->1/2
	],
	LatticeType == "Hypercubic",
	Print @ Show[
		Histogram[
			LatticeEnergies[[All,1,1]],
			{.15}, "PDF",
			ChartStyle->"Pastel"
		],
		ListPlot[
			spectralfunction
		]
	]
]
Print["Integral of spectral function = ", 
	d\[Omega] * Total[spectralfunction[[All,2]]]
];
