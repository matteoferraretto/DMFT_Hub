(* ::Package:: *)

$Path = Join[$Path, {"C:\\Users\\matte\\Desktop\\Mathematica_package\\"}];
(* Import the code *)
<<"DMFT.wl";


(* Import the input file *)
<<"InputFile.wl";


(* Prepare everything *)
<<"Preparation.wl";


(* Start DMFT Loops ... *)
<<"DMFT_Loop.wl";


<<"Post_Processing.wl";
(*KineticEnergyNew[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega], EdMode]*)


(* \[CapitalSigma]0 = Last /@ \[CapitalSigma];

2. * 2. * TMats * Total @ Table[
	LatticeWeights[[k]] * Tr[ (LatticeEnergies[[k]] - \[Mu]) * Total @ Table[
		Inverse[(i\[Omega][[n]] + \[Mu])*IdentityMatrix[Norb] - LatticeEnergies[[k]] - DiagonalMatrix[\[CapitalSigma][[All, n]]]]
		- (LatticeEnergies[[k]] + \[CapitalSigma]0)/(i\[Omega][[n]]^2)
		, {n, CGNMatsubara}]]
	, {k, LatticePoints}] + 
	0.5 * Sum[
		LatticeWeights[[k]] * Tr[(LatticeEnergies[[k]] - \[Mu])]
	, {k, LatticePoints}] - 
	(1./(2.*TMats)) * Sum[
		LatticeWeights[[k]] * Tr[(LatticeEnergies[[k]] - \[Mu]) * (LatticeEnergies[[k]] + \[CapitalSigma]0)]
	, {k, LatticePoints}] *)


(* non interacting kinetic energy *)
N[ 2*2*Integrate[(2/(Pi))*x*Sqrt[1-x^2], {x, -1, 0}] ]


eigs = (Eigensystem[#] &/@ Hsectors);
GF11 = GreenFunctionED[L, f, Norb, {1,1}, {1,1}, {1,1}, Sectors, QnsSectorList, eigs, 0, i\[Omega], EdMode];
GF12 = GreenFunctionED[L, f, Norb, {1,1}, {1,2}, {1,1}, Sectors, QnsSectorList, eigs, 0, i\[Omega], EdMode];

ListPlot[{
	Re[ Gimp[[All,1,2]] ],
	(0.5 + U[[1]]/2.)/(i\[Omega]^2 - (0.5+U[[1]]/2)^2),
	Re[ GF12 ]
}, PlotRange->{{0,500},All}, Joined->True, PlotStyle->{Thick,Dashing[.05],Dashing[.025]}]

ListPlot[{
	Im[ Gimp[[All,1,1]] ],
	Im[ i\[Omega]/(i\[Omega]^2 - (0.5+U[[1]]/2)^2) ],
	Im[ GF11 ]
}, PlotRange->{{0,500},All}, Joined->True, PlotStyle->{Thick,Dashing[.05],Dashing[.025]}]

ListPlot[{
	Re[ \[CapitalSigma][[All,1,2]] ],
	ConstantArray[(U[[1]]/2.), NMatsubara]
}, PlotRange->{{0,500},{0,1.1}}, Joined->True, PlotStyle->{Thick,Dashing[.05]}]


energies = ConstantArray[ ConstantArray[0., {f,f}], {5, Norb, Norb}];
Do[
	(* get the unitary matrix that diagonalizes M: P.M.Pdg = \[CapitalLambda], where \[CapitalLambda] is diagonal *)
	energies[[All, orb, orb]] = (
		(DispersionHypercubicRaman[#, W[[orb]], f, \[Gamma], u] + \[Delta][[orb]]*IdentityMatrix[f] + M[[orb]]) &/@ BrillouinZone[5, 1]
	)
, {orb, Norb}]
energies // Dimensions


hamiltonian[k_, \[CapitalOmega]_] := {
	{-2.*Cos[k], \[CapitalOmega]},
	{\[CapitalOmega], -2.*Cos[k]}
};

\[Epsilon][k_, \[CapitalOmega]_] := Eigenvalues[hamiltonian[k,\[CapitalOmega]]];

Plot[{\[Epsilon][k,2.5][[1]], \[Epsilon][k,2.5][[2]]}, {k,-Pi,Pi}]


Sort @ Eigenvalues[{
{e11, e12, V11, V12},
{e12, e22, V12, V22},
{V11, V12, -\[Lambda], 0},
{V12, V22, 0, \[Lambda]}
} /. {e11 -> 0., e22 -> 0., V11 -> 1., V22 -> 1., V12 -> 0.3, e12 -> 0.3, \[Lambda] -> 2.5}]

Accumulate @ %



G[z_] := Inverse[ (
	z*IdentityMatrix[2] - (-\[Lambda]*PauliMatrix[3]) - {{V11, V12},{V12,V22}} . (Inverse[ {{z-e11, -e12},{-e12,z-e22}}]) . {{V11, V12},{V12,V22}}
) /. {e11 -> 0., e22 -> 0., V11 -> 1., V22 -> 1., V12 -> 0.3, e12 -> 0.3, \[Lambda] -> 2.5}];

ListPlot[{
	Re[ ( G[#] &/@ i\[Omega][[;;CGNMatsubara]] )[[All, 2, 1]] ],
	Re[ Gimp[[All, 2, 1]] ]
}, PlotStyle->{Thick, Dashing[.05]}, Joined->True]


ListPlot[{
	Re[ Gimp[[All, 1,2]] ],
	Re[ LocalG[[All, 1,2]] ]
}, PlotRange->{{0,1000},All}, Joined->True, PlotStyle->{Thick, Dashing[.05]}]


flavdist = Table[
	MomentumDistributedDensityRaman[i, 1, LatticeEnergies, M, \[Mu], \[CapitalSigma], i\[Omega], "Flavor" -> "Real"]
, {i, LatticePoints}];


ListPlot[{
	Re[ flavdist[[All, 1, 1]] ],
	Re[ flavdist[[All, 2, 2]] ]
}, Joined->True, Filling->Axis, PlotRange->{0,1}]

ListPlot[{
	LatticeEnergies[[All, 1, 1, 1, 1]],
	LatticeEnergies[[All, 1, 1, 2, 2]]
}]

2 * FlavorCurrent[W[[1]], \[Gamma], 1, 1, flavdist, LatticeType, LatticeDim, LatticePoints]
2 * FlavorCurrent[W[[1]], \[Gamma], 2, 1, flavdist, LatticeType, LatticeDim, LatticePoints]


ListPlot3D[
	Partition[
		Re[ flavdist[[All, 1, 1]] ],
	31],
	Mesh->None,InterpolationOrder->0,ColorFunction->"SouthwestColors"
]

ListPlot3D[
	Partition[
		LatticeEnergies[[All, 1, 1, 1,2]],
	31],
	Mesh->None,InterpolationOrder->3,ColorFunction->"SouthwestColors"
]


(1./LatticePoints) * Sum[
	Tr[ LatticeEnergies[[i, 1, 1]] . Re[flavdist[[i]]] ]
, {i, LatticePoints}]


Table[
	Mean[ Re[ flavdist[[All, \[Alpha], \[Beta]]] ] ],
{\[Alpha], f}, {\[Beta], f}] // MatrixForm
Total[%]


KineticEnergy[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega], EdMode]


NewBathParameters = BathParameters;


(*\[Omega] = Table[\[Omega]min + n*d\[Omega], {n, 0, NReal}]*)
\[Omega] = Table[-5. + n*0.05, {n, 0, 200}];

eigvecs = Last[ SortBy[Eigensystem[M[[1]]]\[Transpose], First]\[Transpose] ];
Pdg = (Normalize[#] &/@ eigvecs)\[Transpose];
P = ConjugateTranspose[Pdg];

spec = Table[
	- (1./Pi) * Im[
		P . Inverse[ ((\[Omega][[i]] + I*0.05 + \[Mu])*IdentityMatrix[f] - LatticeEnergies[[j, 1, 1]] )] . Pdg
	],
{i, 200}, {j, LatticePoints}];
spec // Dimensions


ListDensityPlot[ 
	spec[[All, All, 1,1]] - spec[[All, All, 2,2]],
	PlotRange -> All,
	ColorFunction -> (RGBColor[If[#>0,#,0], 0, If[#<0, -#, 0]]&),
	ColorFunctionScaling->False
]



max = Max[spec[[All, All, 1,1]]];

Show[{
	ListDensityPlot[
		spec[[All, All, 1,1]]/max,
		PlotRange->All,
		ColorFunction->(Apply[RGBColor, {1,0,0,#}]&),
		ColorFunctionScaling->False,
		DataRange->{{-Pi, Pi}, {\[Omega]min, \[Omega]max}},
		FrameStyle->Directive[Black,16],
		FrameTicks -> {
			{Range[\[Omega]min,\[Omega]max,1.0], None},
			{{{-Pi,"-\[Pi]"}, {-Pi/2, "-\[Pi]/2"}, {0, "0"}, {Pi/2, "\[Pi]/2"}, {Pi, "\[Pi]"}}, None}
		},
		Epilog -> Line[{{-Pi, \[Mu]}, {Pi, \[Mu]}}]
	],
	ListDensityPlot[
		spec[[All, All, 2,2]]/max,
		PlotRange->All,
		ColorFunction->(RGBColor[0,0,1,#]&),
		ColorFunctionScaling->False,
		DataRange->{{-Pi, Pi}, {\[Omega]min, \[Omega]max}},
		FrameStyle->Directive[Black,16],
		FrameTicks -> {
			{Automatic},
			{{{-Pi,"-\[Pi]"}, {-Pi/2, "-\[Pi]/2"}, {0, "0"}, {Pi/2, "\[Pi]/2"}, {Pi, "\[Pi]"}}, None}
		},
		Epilog -> Line[{{-Pi, \[Mu]}, {Pi, \[Mu]}}]
	]
}]


\[Omega] = Table[-7.5 + n*0.05, {n, 0, 300}];

Grealfreq = Mean[Apply[
	GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega] + I*\[Eta]]&,
		{Gs, GsQns}\[Transpose]
	, {1}]]; (* when EdMode == "Raman", the \[Sigma]=1 input is simply ignored, so this is fine! *)
(* G^-1(i\[Omega]) *)
InverseGrealfreq = InverseGreenFunction[Grealfreq, EdMode];
(* G_0^-1(i\[Omega]) *)
InverseG0realfreq = (Weiss/.Thread[symbols -> IndependentParameters])/.{z -> #} &/@ (\[Omega] + I*\[Eta]);
(* \[CapitalSigma](i\[Omega]) *)
\[CapitalSigma]realfreq = InverseG0realfreq - InverseGrealfreq;

spec = MomentumResolvedSpectralFunction[LatticeEnergies[[All, 1,1]], \[Mu], \[CapitalSigma], Range[LatticePoints], \[Omega]+I*\[Eta], EdMode, RamanMatrix -> M[[1]]];
PlotSpectralFunctionRaman[spec, \[Omega], \[Mu], LatticeType, LatticeDim]


<<"Plots.wl";


Show[
ListPlot[
	Transpose[Sort[Eigenvalues[#]] &/@ LatticeEnergies[[All, 1, 1]]]
, DataRange->{-Pi,Pi}],
Plot[
-2*W[[1]]*Cos[\[Gamma][[1,1]]/2]*Cos[k] + M[[1,1,2]]*Sqrt[1 + (2*W[[1]]*Sin[\[Gamma][[1,1]]/2]*Sin[k]/M[[1,1,2]])^2]
, {k,-Pi,Pi}, PlotStyle->Dashed]
]


(* ::Input:: *)
(*{{0.`,2.842170943040401`*^-17},{0.5`,-0.01378565214470845`},{1.`,-0.07149044281346618`},{1.5`,-0.2602015302432129`},{2.`,-0.21392593141245408`},{2.5`,-0.17990361562704466`},{3.`,-0.15444048050075862`},{3.5`,-0.1349090686467734`},{4.`,-0.1195603896542759`},{4.5`,-0.10723293292528117`},{5.`,-0.09714176718383942`},{5.5`,-0.0887440162818923`},{6.`,-0.08165516405378702`}}*)


ListPlot[{
	-Re[ Gimprealfreq[[All,1,2]] ],
	-Re[( Inverse[#] &/@ InverseG0realfreq )[[All, 1,2]] ]
	}, Joined->True, PlotRange->All, PlotStyle->{Thick,Dashed}
]
ListPlot[{
	-Re[InverseG0realfreq[[All,1,1]] ],
	-Re[(Inverse[#] &/@ Gimprealfreq)[[All,1,1]]]
	}, Joined->True, PlotRange->All, PlotStyle->{Thick,Dashed}
]


With[ {\[CapitalOmega] = 0.50, U = 2.},
Show[
	ListPlot[{
		- Im[ Gimprealfreq[[All, 1,2]] ],
		0.5 * Im[ 1./(\[Omega] +I*\[Eta] - \[CapitalOmega] - U/2.) + 1./(\[Omega] +I*\[Eta] + \[CapitalOmega] + U/2.) ]
	}, Joined->True, PlotRange->All, PlotStyle->{Thick,Dashed}]
]
]


\[Omega] = Table[\[Omega]min + d\[Omega]*n, {n, 0, 1000}];
eigs = Eigensystem[#] &/@ Hsectors;
spectralfunctionexact = - (1./Pi) * Im[ GreenFunctionED[L, f, Norb, {1,1}, {1,1}, {1,1}, Sectors, QnsSectorList, eigs, 0, \[Omega] + I*\[Eta], EdMode] ];

ListPlot[{
	{\[Omega], -(1./Pi)*Im[Tr[#] &/@ Gimprealfreq]}\[Transpose],
	{\[Omega], 2*spectralfunctionexact}\[Transpose]
}, Joined->True, PlotRange->All, PlotStyle->{Thick, Dashing[.05]}, AspectRatio->1/2]


Hsectors = SparseArray[#]&/@ With[
{FlatBathParameters = Flatten[Delete[{{0,\[CapitalOmega],0,0,0,0}, \[Delta], {Ug,Ue}, v, v-vex, Jphreshaped, -vex, -mu}, -3]]},
		Sum[
			FlatBathParameters[[i]]*#[[i]],
		{i, 1, Length@FlatBathParameters}] &/@ HlocBlocks
];

QnsSectorList
MatrixForm[#] &/@ Hsectors


eigs = Eigensystem[#] &/@ (Hsectors)/.{\[CapitalOmega]->0.5, v->3., vex->1., Ug->2., Ue->2., mu->(1.+3.-0.5)};

GF = Table[
	GreenFunctionED[L, f, Norb, {1,1}, {\[Sigma], \[Rho]}, {a, b}, Sectors, QnsSectorList, eigs, 0, i\[Omega], EdMode]
, {a, Norb}, {b, Norb}, {\[Sigma],f}, {\[Rho],f}]


ListPlot[
	Re[ GF[[2,2,1,2]] ]
	, Joined->True
]


AtomicGroundState[\[CapitalOmega]_, Ug_, Ue_, V_, Vex_] := Module[
	{FlatBathParameters, Hsectors, eigs, EgsSectorList, Egs, Gs, GsSectorList, GsSectorIndex, GsQns, \[Mu] = Ug/2 + V - Vex/2.},
		FlatBathParameters = Flatten[Delete[{{0,\[CapitalOmega],0,0,0,0}, \[Delta], {Ug,Ue}, V, V-Vex, Jphreshaped, -Vex, -\[Mu]}, -3]];
		Hsectors = Sum[
			FlatBathParameters[[i]]*#[[i]],
		{i, 1, Length@FlatBathParameters}] &/@ HlocBlocks
		(*
		eigs = Eigensystem[#] &/@ Hsectors;
		eigs = SortBy[#, First] &/@ eigs;
		EgsSectorList = eigs[[All, 1]];
		GsSectorList = eigs[[All, -1]];
		Egs = Min[ Flatten[EgsSectorList] ];
		GsSectorIndex = Position[
			EgsSectorList,
			_?((Abs[# - Egs] < DegeneracyThreshold)&)
		];
		(* list of quantum numbers of the degenerate ground states *)
		GsQns = QnsSectorList[[GsSectorIndex[[All, 1]]]];
		If[Ug != Ue, Print["Warning: half filling not guaranteed. "] ];		
		Print["Ground state energy: ", Egs];
		Print["Ground state sector: ", GsQns];
		*)
];

mat = AtomicGroundState[\[CapitalOmega], Ug, Ug, v, Vex];
MatrixForm[#] &/@ mat

FullSimplify[ Eigenvalues[ mat[[5]] ] ]



ListPlot[{
	Re[ \[CapitalSigma][[1,All,1,2]] ],
	Re[ \[CapitalSigma][[2,All,1,2]] ]
},PlotRange->All]

ListPlot[{
	Im[ \[CapitalSigma][[1,All,1,1]] ],
	Im[ \[CapitalSigma][[2,All,1,1]] ]
},PlotRange->All]

Abs[Jse]*0.458 + 3.0*0.456
Abs[Jse]*0.456 + 3.0*0.458

MatrixForm[#] &/@ cdgc
Z


ListPlot[
	Table[
		(Sort[Eigenvalues[#]] &/@ LatticeEnergies[[All,1,1]])[[All,n]],
	{n,4}]
]


Plot[{
	Sort[Eigenvalues[{
		{-2*W[[1]]*Cos[k-Pi/4], M[[1,1,2]]},
		{M[[1,1,2]], -2*W[[1]]*Cos[k+Pi/4]}
	}]],
	Sort[Eigenvalues[{
		{-2*W[[1]]*Cos[k+Pi-Pi/4], M[[1,1,2]]},
		{M[[1,1,2]], -2*W[[1]]*Cos[k+Pi+Pi/4]}
	}]]
	}
	, {k, -Pi/2, Pi/2}
]

Plot[
	Sort[Eigenvalues[{
		{0, 0.25, -W[[1]]*(Exp[I*Pi/4]+Exp[-2*I*k]*Exp[-I*Pi/4]), 0},
		{0.25, 0, 0, -W[[1]]*(Exp[-I*Pi/4]+Exp[-2*I*k]*Exp[I*Pi/4])},
		{-W[[1]]*(Exp[-I*Pi/4]+Exp[2*I*k]*Exp[I*Pi/4]), 0, 0, 0.25},
		{0, -W[[1]]*(Exp[I*Pi/4]+Exp[2*I*k]*Exp[-I*Pi/4]), 0.25, 0}
	}]]
	, {k, -Pi/2, Pi/2}
]

