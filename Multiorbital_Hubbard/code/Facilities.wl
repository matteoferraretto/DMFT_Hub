(* ::Package:: *)

BeginPackage["Facilities`"]

DrawState::usage = "DrawState[L, f, Norb] draws a graphic representation of a Fock state that can be manipulated. Each box can be filled either with 0 (no particles in that slot) or 1 (a particle in that slot). "

EdModeInfo::usage = "EdModeInfo[EdMode] print useful information about EdMode. "

PlotMatsubara::usage = "."

BandPlot2D::usage = "BandPlot2D[Energies, BZ]"

PlotSpectralFunction::usage = "PlotSpectralFunction[spectralfunction]"

PlotSpectralFunctionRaman::usage = "PlotSpectralFunctionRaman[spectralfunction_, \[Omega]_, \[Mu]_, LatticeType_, LatticeDim_]"

WriteOutput::usage = "WriteOutputNew[condition_, OutputDirectory_, label_, data_]"

HNonlocalInfo::usage = "HNonlocalInfo[L, f, Norb, EdMode] print info about the organization of Hamiltonian blocks. "

Wannier1D::usage = "Wannier1D[s, BZ, xlist] computes the lowest-band 1-dimensional Wannier function associated to the lattice potential V(x) = s Er Sin[Pi*x/a]^2, where Er=1 and a=1 are
units of energy and of length respectively. The result is expressed in units of 1/a. BZ is a list of wavevectors in the first Brillouin Zone; while xlist is a list of positions where the
function is computed. An option ''dim'' can be provided: it specifies the dimension of the hamiltonian matrix used in the solution of the central equation; this is 41 by default."


Begin["Private`"];
(*
(* Draw a picture of a state to help the user *)
DrawState[L_, f_, Norb_, j_, \[Sigma]_, orb_] := Module[
	{impuritycolor, color},
	impuritycolor[i_]:=If[i==1,Blue,Black];
	color[i_,index_]:=If[i==j && index==f*(orb-1)+\[Sigma], Red, Black];
	Table[
		Graphics[
			Table[
				{EdgeForm[{Thick,impuritycolor[i]}],color[i,index],Opacity[.3],Rectangle[{1.1*i,0}]},{i,L}
			]
		],
	{index, f*Norb}]
];
DrawState[L_, f_, Norb_] := Module[{},
	Print[Style["Architecture of a state",16]];
	Print[Style["Red:",Red]," (site j, spin \[Sigma], orbital orb)"];
	Print[Style["Blue edge:",Blue]," impurity"];
	Manipulate[
		DrawState[L,f,Norb,j,\[Sigma],orb],
	{j,1,L,1}, {\[Sigma],1,f,1}, {orb,1,Norb,1}]
];

(* info on the EdMode value *)
EdModeInfo[EdMode_] := Which[
	EdMode == "Normal",
	Print["The sectors' quantum numbers are the number of fermions for each flavor and each orbital."],
	EdMode == "Superc",
	Print["The sectors' quantum numbers are the total spin_z operators for each orbital. The bath can exchange pairs with a reservoir, but pairs have an orbital index."],
	EdMode == "InterorbNormal",
	Print["The sectors' quantum numbers are the total number of fermions for each flavor (NOT orbital-wise)."],
	EdMode == "Raman",
	Print["The sectors' quantum numbers are the orbital-wise total number of fermions (NOT flavor-wise)"],
	EdMode == "InterorbSuperc",
	Print["The sectors' quantum numbers are the total spin_z operators (NOT orbital-wise). The bath can exchange pairs with a reservoir, but pairs are inherently interorbital."],
	EdMode == "FullSuperc",
	Print["The sectors' quantum numbers are the total spin_z operators (NOT orbital-wise). The bath can exchange pairs with a reservoir, but pairs are both intraorbital and interorbital."]
];

(* prints useful info about the order of Hamiltonian blocks *)
HNonlocalInfo[L_, f_, Norb_, EdMode_] := Module[{},
	Print["The output blocks have the following order:"]
	Do[
		Print[flag,",  j=",j,",   \[Sigma]=",\[Sigma],",  orb=",orb]
	,{flag,{"Bath","Hopping"}}, {orb,1,Norb}, {\[Sigma],1,f}, {j,2,L}];
	If[
		EdMode == "Superc" || EdMode == "FullSuperc",
		Do[
			Print["Superc,  j=",j,",  orb=",orb]
		,{orb,1,Norb},{j,2,L}];
	];
	If[
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		Do[
			Do[
				If[
					orb2>orb,
					Print["InterorbSuperc,  j=",j,",  orbs=",orb," ",Mod[orb,Norb]+1]
				]
			,{orb2,1,Norb}]
		,{orb,1,Norb},{j,2,L}];
	];
];

(* Show how the Green function is computed within the Nambu formalism *)
GreenFunctionNambuInfo[Norb_] := Module[
	{m = ConstantArray[0, {2Norb, 2Norb}]},
	(* top left element of each block *)
	Do[
		m[[2(orb1-1)+1,2(orb2-1)+1]]=1
	,{orb1,Norb},{orb2,orb1,Norb}];
	(* bottom right element of each block *)
	Do[
		m[[2(orb1-1)+2,2(orb2-1)+2]]=2
	,{orb1,Norb},{orb2,orb1,Norb}];
	(* top right element of each block *)
	Do[
		m[[2(orb1-1)+1,2(orb2-1)+2]]=-1
	,{orb1,Norb},{orb2,orb1,Norb}];
	(* bottom left element of each non-diagonal block *)
	Do[
		m[[2(orb1-1)+2,2(orb2-1)+1]]=-2
	,{orb1,Norb},{orb2,orb1+1,Norb}];
	m//MatrixPlot
];
*)


(* Plot many body function of Matsubara frequencies *)
(*
PlotMatsubara[function_, i\[Omega]_, EdMode_, OptionsPattern[ListPlot]] := Which[
	EdMode == "Normal",
	ListPlot[
		{Abs[i\[Omega]], function}\[Transpose],
		Joined -> True,
		PlotStyle -> {Thick},
		AxesLabel -> {"\!\(\*SubscriptBox[\(\[Omega]\), \(n\)]\)", None}
	],
	EdMode == "Superc",
	ListPlot[{
		Labeled[{Abs[i\[Omega]], function[[All, 1, 1]]}\[Transpose] , "normal"],
		Labeled[{Abs[i\[Omega]], function[[All, 1, 2]]}\[Transpose] , "anomalous"]
		},
		Joined -> True,
		PlotStyle -> {Thick, Thick},
		PlotLayout -> "Row",
		AxesLabel -> {"\!\(\*SubscriptBox[\(\[Omega]\), \(n\)]\)", None}
	]
];
*)

(* plot 2d band structure *)
BandPlot2D[Energies_, BZ_] := Module[
	{energies = Sort[Eigenvalues[#]] &/@ Energies},
		ListPlot3D[
			Table[
				Join[BZ[[#]], {energies[[#, bandindex]]}] &/@ Range[Length[BZ]]
			, {bandindex, Length[energies[[1]]]}],
		(* options *)
			AxesLabel -> {"\!\(\*SubscriptBox[\(k\), \(x\)]\)","\!\(\*SubscriptBox[\(k\), \(y\)]\)","energy"},
			AxesStyle -> Directive[Black, 14],
			PlotStyle -> Opacity[.5]
		]
];

(* plot the spectral function A(\[Omega]) *)
PlotSpectralFunction[spectralfunction_] := Module[{
	rank = ArrayDepth[spectralfunction], \[Omega]min, \[Omega]max, Amax
	},
	Which[
		rank == 2,
		\[Omega]min = spectralfunction[[1,1]];
		\[Omega]max = spectralfunction[[-1,1]];
		Amax = Max[spectralfunction[[All,2]]];,
	(* ----------------------------------- *)
		rank == 3,
		\[Omega]min = Min[spectralfunction[[All,1,1]]];
		\[Omega]max = Max[spectralfunction[[All,-1,1]]];
		Amax = Max[Flatten[spectralfunction[[All,All,2]]]];
	];
	Print @ ListPlot[
		spectralfunction, 
		Joined -> True, 
		Filling -> Axis,
		AspectRatio -> 0.6,
		Axes -> False, Frame -> True,
		FrameStyle -> Directive[Black, 18],
		PlotRange-> {{\[Omega]min - 0.1, \[Omega]max + 0.1}, All},
		FrameLabel -> {"\[HBar]\[Omega]/D", "spectral weight"},
		Epilog -> {Dashing[.05], Line[{{0,0}, {0,Amax}}]}
	];
];

(* Plot flavor resolved spectral function for Raman *)
PlotSpectralFunctionRaman[spectralfunction_, \[Omega]_, \[Mu]_, LatticeType_, LatticeDim_] := Module[
	{\[Omega]min = Min[\[Omega]], \[Omega]max = Max[\[Omega]], f = Length[spectralfunction[[1,1]]], ticks, max, Fermiline, datarange},
	(*If[f > 3, Return["Not supported. Only f=2,3 are available. "] ];*)
	max = Max[spectralfunction[[All, All]]]; (* max value of spin resolved spectral function (normalization factor) *)
	ticks = {
		{Range[\[Omega]min, \[Omega]max, 1.0]/0.25, None},
		Which[
			LatticeDim == 1,
			{{{-Pi,"-\[Pi]"}, {-Pi/2, "-\[Pi]/2"}, {0, "0"}, {Pi/2, "\[Pi]/2"}, {Pi, "\[Pi]"}}, None},
			LatticeDim > 1,
			{None, None}
		]
	};
	datarange = {
		Which[
			LatticeDim == 1, 
			{-Pi, Pi},
			LatticeDim > 1,
			{0, 1}
		],
		{(\[Omega]min + 0.5)/0.25, (\[Omega]max - 0.4)/0.25}
	};
	Fermiline = Which[
		LatticeDim == 1,
		Line[{{-Pi, \[Mu]}, {Pi, \[Mu]}}],
		LatticeDim > 1,
		Line[{{0, \[Mu]}, {1, \[Mu]}}]
	];
	(* plot everything *)
	Show[
		Table[
			ListDensityPlot[
				spectralfunction[[All, All, \[Sigma], \[Sigma]]] / max,
				PlotRange -> All,
				FrameLabel -> {"k","\[Omega]"},
				ColorFunction -> (Apply[RGBColor, Flatten[{ If[f==2 && \[Sigma]==2, {0,0,1}, (*else*) UnitVector[3, \[Sigma]]], #}]]&),
				ColorFunctionScaling -> False,
				DataRange -> datarange,
				FrameStyle -> Directive[Black,20],
				FrameTicks -> ticks, 
				AspectRatio -> 1/GoldenRatio,
				Epilog -> Fermiline
			]
		, {\[Sigma], f}]
	]
];

WriteOutput[condition_, OutputDirectory_, label_, data_] := 
If[
	condition,
	Export[OutputDirectory<>label<>".m", data];
	Print["Data stored on file: "<>label<>".m"];
];

(* compute 1D Wannier function for the potential s Er Sin[Pi*x/a]^2 in units of 1/a *)
Wannier1D[s_, BZ_, xlist_, OptionsPattern[]] := Module[{
	Nx = Length[xlist],
	dx = xlist[[2]]-xlist[[1]],
	a = 1.0,
	Nq = Length[BZ],
	dq = BZ[[2]]-BZ[[1]],
	dim = OptionValue["dim"],
	M, Clist, \[Theta], \[Psi]},
	(* avoid stupid bugs: dim must be odd! *)
	If[EvenQ[dim], dim+=1; Print["warning: dim modified to ",dim]];
	(* Nq x dim x dim *)
	M = SparseArray[{
		Band[{1,1}]->Table[(2j+#*a/Pi)^2+s/2, {j,-(dim-1)/2,(dim-1)/2}],
		Band[{1,2}]->-s/4,
		Band[{2,1}]->-s/4
	}, {dim, dim}] &/@ BZ; 
	(* Nq x dim *)
	Clist = MinimalBy[
		Eigensystem[#, Method->"Banded"]\[Transpose],
		First
	][[1,2]] &/@ M;
	(* Nq *)
	\[Theta] = Arg[Total[#]] &/@ Clist;
	(* vector of dim. Nx *)
	(1.0/Sqrt[2.Pi])*Sqrt[a/(2.0*Pi)]*dq*Total[
		Table[
			Exp[-I*\[Theta][[qindex]]]*Clist[[qindex]] . Table[
				Exp[I*((Pi/a)*2j+BZ[[qindex]])*#]
			, {j, -(dim-1)/2, (dim-1)/2}] 
		,{qindex, Nq}] 
	]&/@ xlist
];
Options[Wannier1D] = {"dim"->41};

(* --- in progress --- *)
(*
BathPlot[L_, BathParameters_] := Module[
	{coordinates, rules, edgeweights},
	coordinates = Join[
		{{0, 0.2}}, (* impurity *)
		{BathParameters[[1,1]], ConstantArray[0, L-1]}\[Transpose]
	];
	edgeweights = BathParameters[[2,1]];
	rules = Thread[ConstantArray[0, L-1] \[UndirectedEdge] Table[n, {n, L-1}]];
	Graph[
		rules,
		VertexCoordinates -> coordinates,
		EdgeStyle -> Thickness[#]&/@Abs[0.25*edgeweights],
		VertexStyle -> Orange
	]
];
*)
Print["Package Facilities` loaded successfully."];

End[];

EndPackage[];
