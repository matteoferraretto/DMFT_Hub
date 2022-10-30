(* ::Package:: *)

BeginPackage["DMFTFatilities`"]

DrawState::usage = "DrawState[L, f, Norb] draws a graphic representation of a Fock state that can be manipulated. Each box can be filled either with 0 (no particles in that slot) or 1 (a particle in that slot). "
EdModeInfo::usage = "EdModeInfo[EdMode] print useful information about EdMode. "
PlotMatsubara::usage = "."
WriteOutput::usage = "WriteOutput[condition, file, label, U, data] "
HNonlocalInfo::usage = "HNonlocalInfo[L, f, Norb, EdMode] print info about the organization of Hamiltonian blocks. "

Begin["`Private`"];

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

(* Plot many body function of Matsubara frequencies *)
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

(* manage output  *)
WriteOutput[condition_, file_, label_, U_, data_] := Module[
	{fout},
	Which[
		label == "BathParameters",
		If[condition,
			Export[file, data, "Table"];
			Print["Converged bath parameters stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label == "SelfEnergy",
		If[condition,
			Export[file, data, "Table"];
			Print["Self Energy stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label == "Error",
		If[condition,
			Export[file, data, "Table"];
			Print["Error list stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label == "SpectralFunction",
		If[condition,
			Export[file, data, "Table"];
			Print["Spectral Function stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label == "z",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, U, " ", DecimalForm@data, "\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["z = ", data, "\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "\[Phi]",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, U, " ", DecimalForm@data, "\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["\[Phi] = ", data, "\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "Ekin",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, U, " ", DecimalForm@data, "\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["\!\(\*SubscriptBox[\(E\), \(kin\)]\) = ", data, "\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "Ds",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, U, " ", DecimalForm@data, "\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["\!\(\*SubscriptBox[\(D\), \(s\)]\) = ", data, "\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label == "Occupancy",
		If[condition,
			fout = OpenAppend[file];(* open output stream on file *)
			WriteString[fout, U, " ", DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]]," ",data[[2]]+2.*data[[3]],"\n"]; (*write the data string*)
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

(* --- in progress --- *)
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

End[];

EndPackage[];
