(* ::Package:: *)

BeginPackage["Parameters`"]

StartingBath::usage = "StartingBath[L, f, Norb, InitializeBathMode, EdMode] returns a list containing the bath parameters to start the DMFT loop.
If EdMode = ''Normal'' then the output has the form {e,V}, where e and V are lists of Norb x (L-1) elements, representing the bath energies and the bath-impurity hybridizations.
If EdMode = ''Superc'' then the output has the form {e,V,\[CapitalDelta]}, where e and V are defined as above, and \[CapitalDelta] is the Norb x Nbath dimensional list of pairs creation (annihilation) amplitudes.
If EdMode = ''InterorbSuperc'' then the output has the form {e,V,\[CapitalDelta],\[CapitalXi]}, where e, V, \[CapitalDelta] are as above, and \[CapitalXi] is the Nbath - dimensional list of interorbital pairs creation (annihilation) amplitudes
InitializeBathMode is a string with the path to the file containing the bath parameters; if it is set to ''Default'', default parameters are dropped."

Symbols::usage = "Symbols[L, f, EdMode] returns a list of symbols representing the independent bath parameters. 
The word ''independent'' means that in case there is some symmetry, for example orbital and spin symmetry, we just extract a representative subset of the bath parameters. "

TakeIndependentParameters::usage = "TakeIndependentParameters[L, f, Norb, \[Sigma], orb, BathParameters, EdMode] returns a flat list of independent bath parameters depending on EdMode.
The word ''independent'' means that in case there is some symmetry, for example orbital and spin symmetry, we just extract a representative subset of the bath parameters with labels 
\[Sigma] and orb, namely e_\[Sigma],orb ; V_\[Sigma],orb, etc. This is useful in two cases: to compute the numerical Weiss field from the symbolic expression, and to perform the self consistency minimization
with the smallest possible number of variables. "

ReshapeBathParameters::usage = "ReshapeBathParameters[L, f, Norb, IndependentParameters, OrbitalSymmetry, EdMode] takes a list of independent bath parameters and returns a reshaped version of it which is 
consistent with the conventionally chosen shape. "

Begin["Private`"]

Print["Package Parameters` loaded successfully."];

(* Initialize starting bath *)
StartingBath[L_, f_, Norb_, InitializeBathMode_, EdMode_, OptionsPattern[]] := Module[
	{e, V, \[CapitalDelta], \[CapitalXi], Nbath},
	Nbath = L - 1;
	Which[
		EdMode == "Normal" || EdMode == "InterorbNormal",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb],	
		(*else*)
			{e,V} = Import[InitializeBathMode,"Table"];
		];
		Return[{e, V}],
(* ---------------------------------------------- *)
		EdMode == "Raman",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[
					ConstantArray[OptionValue[\[CapitalOmega]0], {f, f}] + 
					(-OptionValue[\[CapitalOmega]0]-(Nbath-1)/2.+k)*IdentityMatrix[f]
				, {k, 0, Nbath-1}]
			, Norb];
			V = ConstantArray[
				Table[
					ConstantArray[1., {f, f}]
				, {k, Nbath}]
			, Norb],
		(*else*)
			{e, V} = Import[InitializeBathMode, "Table"];
		];
		Return[{e, V}],
(* ---------------------------------------------- *)
		EdMode == "Superc",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k, {k, 0, Nbath-1}],
			f*Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb];
			\[CapitalDelta] = OptionValue[\[CapitalDelta]0] *ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb],
		(*else*)
			{e,V,\[CapitalDelta]} = Import[InitializeBathMode,"Table"];
		];
		Return[{e, V, \[CapitalDelta]}],
(* ---------------------------------------------- *)
		EdMode == "InterorbSuperc",
		If[
			InitializeBathMode=="Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb];
			V = ConstantArray[
				Table[1., {k,1,Nbath}],
			f*Norb];
	        \[CapitalXi] = OptionValue[\[CapitalXi]0] * Table[1., {k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalXi]} = Import[InitializeBathMode,"Table"];
		];
	   Return[{e,V,\[CapitalXi]}],
(* ---------------------------------------------- *)
		EdMode == "FullSuperc",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb];
			\[CapitalDelta] = OptionValue[\[CapitalDelta]0] * ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb];
	        \[CapitalXi] = OptionValue[\[CapitalXi]0] * Table[1.,{k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalDelta],\[CapitalXi]} = Import[InitializeBathMode,"Table"];
		];
	   Return[{e,V,\[CapitalDelta],\[CapitalXi]}];
	]
];
Options[StartingBath] = {\[CapitalDelta]0 -> 1., \[CapitalXi]0 -> 1., \[CapitalOmega]0 -> 1.};

(* generate a list of "independent" symbols *)
Symbols[L_, f_, Norb_, EdMode_] := Which[
	EdMode == "Normal", 
	Join[
		Table[Symbol["e"<>ToString[i]], {i, L-1}],
		Table[Symbol["V"<>ToString[i]], {i, L-1}]
	],
(* ------------------------------------------------ *)
	EdMode == "Superc",
	Join[
		Table[Symbol["e"<>ToString[i]], {i, L-1}],
		Table[Symbol["V"<>ToString[i]], {i, L-1}],
		Table[Symbol["\[CapitalDelta]"<>ToString[i]], {i, L-1}]
	],
(* ------------------------------------------------ *)
	EdMode == "Raman",
	Join[
		Flatten[
			Table[
				Symbol["e"<>ToString[i]<>ToString[n]<>ToString[m]]
			, {i, 1, L-1}, {m, 1, f}, {n, m, f}]
		, 3],
		Flatten[
			Table[
				Symbol["V"<>ToString[i]<>ToString[n]<>ToString[m]]
			, {i, 1, L-1}, {m, 1, f}, {n, m, f}]
		, 3]
	],
(* ---------------------------------------------- *)
	EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
	Join[
		(* e11, e12, e13 ... -> orb = 1 site i = 1,2,3... *)
		Flatten @ Table[
			Symbol["e"<>ToString[orb]<>ToString[i]]
		, {orb, Norb}, {i, L-1}],
		(* V11, V12, V13 ... -> orb = 1 site i = 1,2,3... *)
		Flatten @ Table[
			Symbol["V"<>ToString[orb]<>ToString[i]]
		, {orb, Norb}, {i, L-1}],
		(* \[CapitalDelta]11, \[CapitalDelta]12, \[CapitalDelta]13 ... -> orb = 1 site i = 1,2,3... *)
		Flatten @ Table[
			Symbol["\[CapitalDelta]"<>ToString[orb]<>ToString[i]]
		, {orb, Norb}, {i, L-1}],
		(* \[CapitalXi]1, \[CapitalXi]2, \[CapitalXi]3 ... where 1,2,3... are site indexes *)
		Table[
			Symbol["\[CapitalXi]"<>ToString[i]]
		, {i, L-1}]
	]
];

(* extract independent bath parameters depending on EdMode *)
TakeIndependentParameters[L_, f_, Norb_, \[Sigma]_, orb_, BathParameters_, EdMode_] :=
	Which[
		EdMode == "Normal",
		Join[
			BathParameters[[1]][[f*(orb-1)+\[Sigma]]], (* e1, e2, e3, ... *)
			BathParameters[[2]][[f*(orb-1)+\[Sigma]]] (* V1, V2, V3, ... *)
		],
	(* ---------------------------------------- *)
		EdMode == "Superc",
		Join[
			BathParameters[[1]][[f*(orb-1)+\[Sigma]]], (* e1, e2, e3, ... *)
			BathParameters[[2]][[f*(orb-1)+\[Sigma]]], (* V1, V2, V3, ... *)
			BathParameters[[3]][[orb]] (* \[CapitalDelta]1, \[CapitalDelta]2, \[CapitalDelta]3, ... *)
		],
	(* --------------------------------------- *)
		EdMode == "Raman",
		Flatten @ Join[
			Table[
				Join @@ Pick[
					BathParameters[[1]][[orb]][[i]]
				, UpperTriangularize[ConstantArray[1, {f, f}]], 1] (* upper triangular part of e_orb,i *)
			, {i, L-1}],
			Table[
				Join @@ Pick[
					BathParameters[[2]][[orb]][[i]]
				, UpperTriangularize[ConstantArray[1, {f, f}]], 1] (* upper triangular part of V_orb,i *)
			, {i, L-1}]
		],
	(* --------------------------------------- *)
		EdMode == "InterorbNormal", 
		Flatten[BathParameters] (* <- to be fixed *),
	(* --------------------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		Join[
			Flatten @ Table[
				BathParameters[[1]][[f*(o-1)+\[Sigma]]]
			, {o, Norb}], (* e11, e12, ..., e21, e22, ... *)
			Flatten @ Table[
				BathParameters[[2]][[f*(o-1)+\[Sigma]]]
			, {o, Norb}], (* V11, V12, ..., V21, V22, ... *)
			Flatten @ BathParameters[[3]], (* \[CapitalDelta]11, \[CapitalDelta]12, ..., \[CapitalDelta]21, \[CapitalDelta]22, ... *)
			BathParameters[[4]] (* \[CapitalXi]1, \[CapitalXi]2, ... *)
		]
	];

(* correctly reshape the flat list of independent bath parameters *)
ReshapeBathParameters[L_, f_, Norb_, IndependentParameters_, OrbitalSymmetry_, EdMode_] := 
	Which[
		EdMode == "Normal" && OrbitalSymmetry,
		{ConstantArray[Take[IndependentParameters, L-1], f*Norb], (* e *)
		ConstantArray[Take[IndependentParameters, {L, 2(L-1)}], f*Norb]}, (* V *)
	(* ----------------------------------------- *)
		EdMode == "Normal" && !OrbitalSymmetry,
		{Join[Flatten[
			Table[
				ConstantArray[IndependentParameters[[orb]][[;;L-1]], f],
			{orb, Norb}]
		, 1]],
		Join[Flatten[
			Table[
				ConstantArray[IndependentParameters[[orb]][[L;;]], f],
			{orb, Norb}]
		, 1]]},
	(* ----------------------------------------- *)
		EdMode == "Superc" && OrbitalSymmetry,
		{ConstantArray[Take[IndependentParameters, L-1], f*Norb], (* e *)
		ConstantArray[Take[IndependentParameters, {L, 2(L-1)}], f*Norb], (* V *)
		ConstantArray[Take[IndependentParameters, {2L-1, 3(L-1)}], Norb]}, (* \[CapitalDelta] *)
	(* ----------------------------------------- *)
		EdMode == "Superc" && !OrbitalSymmetry,
		{Join[Flatten[
			Table[
				ConstantArray[IndependentParameters[[orb]][[;;L-1]], f],
			{orb, Norb}]
		, 1]],
		Join[Flatten[
			Table[
				ConstantArray[IndependentParameters[[orb]][[L;;2(L-1)]], f],
			{orb, Norb}]
		, 1]],
		Table[
			IndependentParameters[[orb]][[2L-1;;]],
		{orb, Norb}]
		},
	(* ----------------------------------------- *)
		EdMode == "Raman" && OrbitalSymmetry,
		Return[0],
	(* ----------------------------------------- *)
		EdMode == "Raman" && !OrbitalSymmetry,
		Return[0],	
	(* ----------------------------------------- *)
		EdMode == "InterorbNormal",
		{Partition[Take[IndependentParameters, (L-1)*f*Norb], L-1], (* e *)
		Partition[Take[IndependentParameters, {1+(L-1)*f*Norb, 2*(L-1)*f*Norb}], L-1]}, (* V *)
	(* ----------------------------------------- *)
		EdMode == "InterorbSuperc",
		{Partition[Take[IndependentParameters, (L-1)*f*Norb], L-1], (* e *)
		Partition[Take[IndependentParameters, {1+(L-1)*f*Norb, 2*(L-1)*f*Norb}], L-1], (* V *)
		Take[IndependentParameters, {1+2*(L-1)*f*Norb, 2*(L-1)*f*Norb + L-1}]}, (* \[CapitalXi] *)
	(* ----------------------------------------- *)
		EdMode == "FullSuperc",
		{
		Join @@ (ConstantArray[#, f] &/@
			Partition[IndependentParameters[[;;(L-1)*Norb]], L-1]
		), (* {{e11, e12, ...}, {e11, e12, ...}, ... (f times), {e21, e22, ...}}, {e21, e22, ...}}, ... (f times), ... *)
		Join @@ (ConstantArray[#, f] &/@
			Partition[IndependentParameters[[(L-1)*Norb+1 ;; 2*(L-1)*Norb]], L-1]
		), (* {{V11, V12, ...}, {V11, V12, ...}, ... (f times), {V21, V22, ...}}, {V21, V22, ...}}, ... (f times), ... *)
		Partition[IndependentParameters[[2*(L-1)*Norb+1 ;; 3*(L-1)*Norb]], L-1], (* {{\[CapitalDelta]11, \[CapitalDelta]12, ...}, {\[CapitalDelta]21, \[CapitalDelta]22, ...}, ... (Norb times)} *)
		IndependentParameters[[3*(L-1)*Norb+1 ;;]] (* {\[CapitalXi]1, \[CapitalXi]2, ...} *)
		}
	];

End[]

EndPackage[]
