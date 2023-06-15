(* ::Package:: *)

BeginPackage["Facilities`"]

DrawState::usage = "DrawState[L, f, Norb] draws a graphic representation of a Fock state that can be manipulated. Each box can be filled either with 0 (no particles in that slot) or 1 (a particle in that slot). "

EdModeInfo::usage = "EdModeInfo[EdMode] print useful information about EdMode. "

PlotMatsubara::usage = "."

BandPlot2D::usage = "BandPlot2D[Energies, BZ]"

PlotSpectralFunctionRaman::usage = "PlotSpectralFunctionRaman[spectralfunction_, \[Omega]_, \[Mu]_, LatticeType_, LatticeDim_]"

WriteOutput::usage = "WriteOutputNew[condition_, OutputDirectory_, label_, data_]"

HNonlocalInfo::usage = "HNonlocalInfo[L, f, Norb, EdMode] print info about the organization of Hamiltonian blocks. "

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

(* Plot flavor resolved spectral function for Raman *)
PlotSpectralFunctionRaman[spectralfunction_, \[Omega]_, \[Mu]_, LatticeType_, LatticeDim_] := Module[
	{\[Omega]min = Min[\[Omega]], \[Omega]max = Max[\[Omega]], f = Length[spectralfunction[[1,1]]], ticks, max, Fermiline, datarange},
	If[f > 3, Return["Not supported. Only f=2,3 are available. "] ];
	max = Max[spectralfunction[[All, All]]]; (* max value of spin resolved spectral function (normalization factor) *)
	ticks = {
		{Range[\[Omega]min, \[Omega]max, 1.0], None},
		Which[
			LatticeDim == 1,
			{{{-Pi,"-\[Pi]"}, {-Pi/2, "-\[Pi]/2"}, {0, "0"}, {Pi/2, "\[Pi]/2"}, {Pi, "\[Pi]"}}, None}
		]
	};
	datarange = {
		Which[
			LatticeDim == 1, 
			{-Pi, Pi}
		],
		{\[Omega]min, \[Omega]max}
	};
	Fermiline = Which[
		LatticeDim == 1,
		Line[{{-Pi, \[Mu]}, {Pi, \[Mu]}}]
	];
	(* plot everything *)
	Show[
		Table[
			ListDensityPlot[
				spectralfunction[[All, All, \[Sigma], \[Sigma]]] / max,
				PlotRange -> All,
				FrameLabel -> {"k","\[Omega]"},
				FrameStyle -> Directive[Black, 18],
				ColorFunction -> (Apply[RGBColor, Flatten[{ If[f==2 && \[Sigma]==2, {0,0,1}, (*else*) UnitVector[3, \[Sigma]]], #}]]&),
				ColorFunctionScaling -> False,
				DataRange -> datarange,
				FrameStyle->Directive[Black,16],
				FrameTicks -> ticks, 
				Epilog -> Fermiline
			]
		, {\[Sigma], f}]
	]
];

(* manage output  
WriteOutput[condition_, file_, label_, U_, data_] := Module[
	{fout},
	Which[
		label == "BathParameters",
		If[condition,
			Export[file, data];
			(*Print["Converged bath parameters stored on file."]*)
		],
	(* ---------------------------------------- *)
		label == "SelfEnergy",
		If[condition,
			Export[file, data];
			Print["Self Energy stored on file."]
		],
	(* ---------------------------------------- *)
		label == "Error",
		If[condition,
			Export[file, data];
			(*Print["Error list stored on file."]*)
		],
	(* ---------------------------------------- *)
		label == "SpectralFunction",
		If[condition,
			Export[file, data];
			Print["Spectral Function stored on file."]
		],
	(* ---------------------------------------- *)
		label == "z",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, data, "\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["z saved to file."];
		Print["z = ", data, "\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "\[Phi]",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, data, "\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["\[Phi] = ", data, "\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "Ekin",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, data, "\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["\!\(\*SubscriptBox[\(E\), \(kin\)]\) = ", data, "\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "Ds",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, DecimalForm@data, "\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["\!\(\*SubscriptBox[\(D\), \(s\)]\) = ", data, "\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "Occupancy",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]]," ",data[[2]]+2.*data[[3]],"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["U","\t\t\t","Double Occ.","\t\t\t","Single Occ.","\t\t\t","Empty Occ.","\t\t\t","Density"];
		Print[U,"\t\t\t",data[[1]],"\t\t\t ",data[[2]],"\t\t\t\t",data[[3]],"\t\t\t\t",data[[2]]+2.*data[[3]],"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "Density",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]],"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["U","\t\t\t\t","<n>","\t\t\t\t","<\!\(\*SuperscriptBox[\(n\), \(2\)]\)>","\t\t\t\t","\[CapitalDelta]n"];
		Print[U,"\t\t\t",data[[1]],"\t\t\t ",data[[2]],"\t\t\t\t",data[[3]],"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "Spin",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]],"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["U","\t\t\t\t","<\!\(\*SubscriptBox[\(s\), \(z\)]\)>","\t\t\t\t","<\!\(\*SuperscriptBox[SubscriptBox[\(s\), \(z\)], \(2\)]\)>","\t\t\t","\!\(\*SubscriptBox[\(\[CapitalDelta]s\), \(z\)]\)"];
		Print[U,"\t\t\t",data[[1]],"\t\t\t ",data[[2]],"\t\t\t\t",data[[3]],"\n"](* print data on screen *)
	]
];
*)

WriteOutput[condition_, OutputDirectory_, label_, data_] := 
If[
	condition,
	Export[OutputDirectory<>label<>".m", data];
	Print["Data stored on file: "<>label<>".m"];
];

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
