(* ::Package:: *)

BeginPackage["Parameters`", {"MyLinearAlgebra`"}]

StartingBath::usage = "StartingBath[L, f, Norb, \[Delta], InitializeBathMode, EdMode] returns a list containing the bath parameters to start the DMFT loop.
If EdMode = ''Normal'' then the output has the form {e,V}, where e and V are lists of Norb x (L-1) elements, representing the bath energies and the bath-impurity hybridizations.
If EdMode = ''Superc'' then the output has the form {e,V,\[CapitalDelta]}, where e and V are defined as above, and \[CapitalDelta] is the Norb x Nbath dimensional list of pairs creation (annihilation) amplitudes.
If EdMode = ''Raman'' then the output has the form {e, V}, where both e and V are lists with dimension (Norb x Nbath), where each element is a fxf matrix. 
The diagonal components of e[[orb, i]] are on-site flavor-resolved energies for the bath sites; while the off-diagonal components are hybridizations between different flavors 
of the same bath site. The diagonal components of V[[i]] are flavor-diagonal tunnelings between the impurity and the i-th bath site; while the off-diagonal components are the 
flavor-swapping tunnelings between the impurity and the i-th site.
If EdMode = ''InterorbSuperc'' then the output has the form {e,V,\[CapitalDelta],\[CapitalXi]}, where e, V, \[CapitalDelta] are as above, and \[CapitalXi] is the Nbath - dimensional list of interorbital pairs creation (annihilation) amplitudes
InitializeBathMode is a string with the path to the file containing the bath parameters; if it is set to ''Default'', default parameters are dropped.
The input \[Delta] is the list of crystal field splittings (with Norb elements) and it's used to set the starting bath energy levels around the corresponding non-interacting energies."

LocalParameters::usage = "LocalParameters[M, \[Delta], U, Ust, Usec, Jph, Jse, \[Mu], shift, SublatticesQ, hAFM, V] returns a list of local coupling constants in suitable order and shape to be
correctly combined with the Hamiltonian blocks. The meaning of all the symbols is specified in the input file. "

Symbols::usage = "Symbols[L, f, Norb, EdMode] returns a list of symbols representing the independent bath parameters. 
The word ''independent'' means that in case there is some symmetry, for example orbital and spin symmetry, we just extract a representative subset of the bath parameters. "

SymmetrizeSymbols::usage = "SymmetrizeSymbols[symbols_, EdMode_, OrbitalSymmetry_, ParticleHoleSymmetry_]"

IndependentSymbols::usage = ""

IndependentSymbolsIndexes::usage = ""

TakeIndependentParameters::usage = "TakeIndependentParameters[L, f, Norb, \[Sigma], orb, BathParameters, EdMode] returns a flat list of independent bath parameters depending on EdMode.
The word ''independent'' means that in case there is some symmetry, for example orbital and spin symmetry, we just extract a representative subset of the bath parameters with labels 
\[Sigma] and orb, namely e_\[Sigma],orb ; V_\[Sigma],orb, etc. This is useful in two cases: to compute the numerical Weiss field from the symbolic expression, and to perform the self consistency minimization
with the smallest possible number of variables. "

ReshapeBathParameters::usage = "ReshapeBathParameters[L, f, Norb, IndependentParameters, OrbitalSymmetry, EdMode] takes a list of independent bath parameters and returns a reshaped version of it which is 
consistent with the conventionally chosen shape. "

Begin["Private`"]

Print["Package Parameters` loaded successfully."];

(* Initialize starting bath *)
(*StartingBathOld[L_, f_, Norb_, \[Delta]_, InitializeBathMode_, EdMode_, SublatticesQ_, OptionsPattern[]] := Module[
	{e, V, \[CapitalDelta], \[CapitalXi], Nbath},
	Nbath = L - 1;
	Which[
		EdMode == "Normal" || EdMode == "InterorbNormal",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb] + 
			Join @@ Table[
				ConstantArray[\[Delta][[orb]], f]
			, {orb, Norb}];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb] * OptionValue[V0],	
		(*else*)
			{e,V} = Import[InitializeBathMode];
		];
		Return[{e, V}],
(* ---------------------------------------------- *)
		EdMode == "Raman",
		If[InitializeBathMode == "Default",
			e = ConstantArray[
				Table[
					ConstantArray[OptionValue[\[CapitalOmega]0], {f, f}] + 
					(-OptionValue[\[CapitalOmega]0]-(Nbath-1)/2.+k)*IdentityMatrix[f]
				, {k, 0, Nbath-1}]
			, Norb];
			V = ConstantArray[
				Table[
					DiagonalMatrix[ ConstantArray[1.0 - OptionValue[\[CapitalOmega]0], f] ] + 
					ConstantArray[OptionValue[\[CapitalOmega]0], {f, f}]
				, {k, Nbath}]
			, Norb] * OptionValue[V0];
			(* for sublattice calculations duplicate the list *)
			If[SublatticesQ, 
				Return[{{e,V}, {e,V}}],
				Return[{e, V}]
			];,
		(* else, if you import from file *)
			Return[ Import[InitializeBathMode] ];
		],
(* ---------------------------------------------- *)
		EdMode == "Superc",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k, {k, 0, Nbath-1}],
			f*Norb] + 
			Join @@ Table[
				ConstantArray[\[Delta][[orb]], f]
			, {orb, Norb}];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb] * OptionValue[V0];
			\[CapitalDelta] = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb] * OptionValue[\[CapitalDelta]0],
		(*else*)
			{e,V,\[CapitalDelta]} = Import[InitializeBathMode];
		];
		Return[{e, V, \[CapitalDelta]}],
(* ---------------------------------------------- *)
		EdMode == "InterorbSuperc",
		If[
			InitializeBathMode=="Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb] + 
			Join @@ Table[
				ConstantArray[\[Delta][[orb]], f]
			, {orb, Norb}];
			V = ConstantArray[
				Table[1., {k,1,Nbath}],
			f*Norb] * OptionValue[V0];
	        \[CapitalXi] = OptionValue[\[CapitalXi]0] * Table[1., {k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalXi]} = Import[InitializeBathMode];
		];
	   Return[{e,V,\[CapitalXi]}],
(* ---------------------------------------------- *)
		EdMode == "FullSuperc",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb] + 
			Join @@ Table[
				ConstantArray[\[Delta][[orb]], f]
			, {orb, Norb}];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb] * OptionValue[V0];
			\[CapitalDelta] = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb] * OptionValue[\[CapitalDelta]0];
	        \[CapitalXi] = OptionValue[\[CapitalXi]0] * Table[1.,{k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalDelta],\[CapitalXi]} = Import[InitializeBathMode];
		];
	    Return[{e,V,\[CapitalDelta],\[CapitalXi]}];
	]
];
Options[StartingBathOld] = {V0 -> 1., \[CapitalDelta]0 -> 1., \[CapitalXi]0 -> 1., \[CapitalOmega]0 -> 1.};*)

StartingBath[L_, f_, Norb_, \[Delta]_, InitializeBathMode_, EdMode_, SublatticesQ_, OptionsPattern[]] := Module[
	{e, V, \[CapitalDelta], \[CapitalXi], Nbath},
	Nbath = L - 1;
	Which[
		EdMode == "Normal" || EdMode == "InterorbNormal",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			Norb] + \[Delta];
			V = ConstantArray[OptionValue[V0], {Norb, Nbath}],
		(*else*)
			{e,V} = Import[InitializeBathMode];
		];
		Return[{e, V}],
(* ---------------------------------------------- *)
		EdMode == "Raman",
		If[InitializeBathMode == "Default",
			e = ConstantArray[
				Table[
					ConstantArray[OptionValue[\[CapitalOmega]0], {f, f}] + 
					(-OptionValue[\[CapitalOmega]0]-(Nbath-1)/2.+k)*IdentityMatrix[f]
				, {k, 0, Nbath-1}]
			, Norb];
			V = ConstantArray[
				Table[
					DiagonalMatrix[ ConstantArray[1.0 - OptionValue[\[CapitalOmega]0], f] ] + 
					ConstantArray[OptionValue[\[CapitalOmega]0], {f, f}]
				, {k, Nbath}]
			, Norb] * OptionValue[V0];
			(* for sublattice calculations duplicate the list *)
			If[SublatticesQ, 
				Return[{{e,V}, {e,V}}],
				Return[{e, V}]
			];,
		(* else, if you import from file *)
			Return[ Import[InitializeBathMode] ];
		],
(* ---------------------------------------------- *)
		EdMode == "Magnetic",
		If[InitializeBathMode == "Default",
			e = ConstantArray[
				Table[
					ConstantArray[-(Nbath-1)/2. + k, f]
				, {k, 0, Nbath-1}]
			, Norb];
			V = ConstantArray[1.0, {Norb, Nbath, f}] * OptionValue[V0];
			(* for sublattice calculations duplicate the list *)
			If[SublatticesQ, 
				Return[{{e,V}, {e,V}}],
				Return[{e, V}]
			];,
		(* else, if you import from file *)
			Return[ Import[InitializeBathMode] ];
		],
(* ---------------------------------------------- *)
		EdMode == "Superc",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k, {k, 0, Nbath-1}],
			Norb] + \[Delta];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb] * OptionValue[V0];
			\[CapitalDelta] = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb] * OptionValue[\[CapitalDelta]0],
		(*else*)
			{e,V,\[CapitalDelta]} = Import[InitializeBathMode];
		];
		Return[{e, V, \[CapitalDelta]}],
(* ---------------------------------------------- *)
		EdMode == "InterorbSuperc",
		If[
			InitializeBathMode=="Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb] + 
			Join @@ Table[
				ConstantArray[\[Delta][[orb]], f]
			, {orb, Norb}];
			V = ConstantArray[
				Table[1., {k,1,Nbath}],
			f*Norb] * OptionValue[V0];
	        \[CapitalXi] = OptionValue[\[CapitalXi]0] * Table[1., {k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalXi]} = Import[InitializeBathMode];
		];
	   Return[{e,V,\[CapitalXi]}],
(* ---------------------------------------------- *)
		EdMode == "FullSuperc",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb] + 
			Join @@ Table[
				ConstantArray[\[Delta][[orb]], f]
			, {orb, Norb}];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb] * OptionValue[V0];
			\[CapitalDelta] = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb] * OptionValue[\[CapitalDelta]0];
	        \[CapitalXi] = OptionValue[\[CapitalXi]0] * Table[1.,{k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalDelta],\[CapitalXi]} = Import[InitializeBathMode];
		];
	    Return[{e,V,\[CapitalDelta],\[CapitalXi]}];
	]
];
Options[StartingBath] = {V0 -> 1., \[CapitalDelta]0 -> 1., \[CapitalXi]0 -> 1., \[CapitalOmega]0 -> 0.};

(* properly organize local coupling constants *) 
(*LocalParameters[M_, \[Delta]_, U_, Ust_, Usec_, Jph_, Jse_, \[Mu]_, shift_, SublatticesQ_, hAFM_, V_] := Module[
	{Jphreshaped, Mreshaped, MreshapedA, MreshapedB, InteractionParameters, Norb = Length[M]},
	(* just pick upper triangular part of the pair hopping matrix *)
	Jphreshaped = ArrayFromSymmetricMatrix[Jph];
	(* reshape Raman field to match number of corresponding hamiltonian blocks *)
	Mreshaped = Flatten[ArrayFromSymmetricMatrix[#] &/@ M];
	(* get list of interaction parameters in correct order *)
	InteractionParameters = {Mreshaped, \[Delta], U, Ust, Usec, Jphreshaped, Jse, - \[Mu] + shift};
	(* if there are sublattices, duplicate the bath and interaction parameters *)
	If[SublatticesQ,
		InteractionParameters = {InteractionParameters, InteractionParameters};
		If[hAFM != 0 || V != 0,
			MreshapedA = Flatten[ArrayFromSymmetricMatrix[#] &/@ (M + ConstantArray[DiagonalMatrix[{hAFM, -hAFM}], Norb])];
			MreshapedB = Flatten[ArrayFromSymmetricMatrix[#] &/@ (M - ConstantArray[DiagonalMatrix[{hAFM, -hAFM}], Norb])];
			InteractionParameters = {
				{MreshapedA, \[Delta], U, Ust, Usec, Jphreshaped, Jse, - \[Mu] + V + shift},
				{MreshapedB, \[Delta], U, Ust, Usec, Jphreshaped, Jse, - \[Mu] - V + shift}
			};
		]
	];
	InteractionParameters
];*)

LocalParameters[M_, \[Delta]_, U_, Ust_, Usec_, Jph_, Jse_, \[CapitalDelta]ext_, \[Mu]_, shift_, SublatticesQ_, hAFM_, V_, EdMode_] := Module[
	{Mreshaped, MreshapedA, MreshapedB, InteractionParameters, Norb = Length[M]},
	(* reshape Raman field to match number of corresponding hamiltonian blocks *)
	(*Mreshaped = Flatten[ArrayFromSymmetricMatrix[#] &/@ M];*)
	(* get list of interaction parameters in correct order *)
	Which[
		EdMode == "Normal",
		InteractionParameters = {U, Ust, Usec, \[Mu] - shift, \[Delta]},
		EdMode == "Superc",
		InteractionParameters = {U, Ust, Usec, Jph, \[CapitalDelta]ext, \[Mu] - shift, \[Delta]},
		EdMode == "Raman",
		InteractionParameters = {M, U, Ust, Usec, Jse, \[Mu] - shift, \[Delta]} ,
		EdMode == "Magnetic",
		InteractionParameters = {Diagonal[#] &/@ M, U, Ust, Usec, \[Mu] - shift, \[Delta]}
	]
	(* if there are sublattices, duplicate the bath and interaction parameters *)
	If[SublatticesQ,
		InteractionParameters = {InteractionParameters, InteractionParameters};
		If[hAFM != 0 || V != 0,
			MreshapedA = Flatten[ArrayFromSymmetricMatrix[#] &/@ (M + ConstantArray[DiagonalMatrix[{hAFM, -hAFM}], Norb])];
			MreshapedB = Flatten[ArrayFromSymmetricMatrix[#] &/@ (M - ConstantArray[DiagonalMatrix[{hAFM, -hAFM}], Norb])];
			InteractionParameters = {
				{MreshapedA, \[Delta], U, Ust, Usec, Jph, Jse, - \[Mu] + V + shift},
				{MreshapedB, \[Delta], U, Ust, Usec, Jph, Jse, - \[Mu] - V + shift}
			};
		]
	];
	InteractionParameters
];

(* generate a list of "independent" symbols *)
(*Symbols[L_, f_, Norb_, EdMode_] := Which[
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
				Symbol["e"<>ToString[i]<>ToString[m]<>ToString[n]]
			, {i, 1, L-1}, {m, 1, f}, {n, m, f}]
		, 3],
		Flatten[
			Table[
				Symbol["V"<>ToString[i]<>ToString[m]<>ToString[n]]
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
];*)
(* new version *)
Symbols[L_, f_, Norb_, EdMode_, filename_, LoadSymbolsQ_] := Module[
	{e, V, \[CapitalDelta], symbols},
	If[LoadSymbolsQ,
		symbols = Import[filename];,
	(* else *)
		Which[
			EdMode == "Normal",
			e = Table[
				Symbol[ "e"<>"orb"<>ToString[orb]<>"l"<>ToString[l] ]
			, {orb, Norb}, {l, L-1}];
			V = Table[
				Symbol[ "V"<>"orb"<>ToString[orb]<>"l"<>ToString[l] ]
			, {orb, Norb}, {l, L-1}];
			symbols = {e, V};,
		(* ---------------------------------------------------- *)
			EdMode == "Superc",
			e = Table[
				Symbol[ "e"<>"orb"<>ToString[orb]<>"l"<>ToString[l] ]
			, {orb, Norb}, {l, L-1}];
			V = Table[
				Symbol[ "V"<>"orb"<>ToString[orb]<>"l"<>ToString[l] ]
			, {orb, Norb}, {l, L-1}];
			\[CapitalDelta] = Table[
				Symbol[ "\[CapitalDelta]"<>"orb"<>ToString[orb]<>"l"<>ToString[l] ]
			, {orb, Norb}, {l, L-1}];
			symbols = {e, V, \[CapitalDelta]};,
		(* ---------------------------------------------------- *)
			EdMode == "Raman",
			e = Table[
				If[\[Sigma] <= \[Rho],
					Symbol["e"<>"orb"<>ToString[orb]<>"l"<>ToString[l]<>"\[Sigma]"<>ToString[\[Sigma]]<>"\[Rho]"<>ToString[\[Rho]]],
				(* else \[Sigma] < \[Rho] *)
					Symbol["e"<>"orb"<>ToString[orb]<>"l"<>ToString[l]<>"\[Sigma]"<>ToString[\[Rho]]<>"\[Rho]"<>ToString[\[Sigma]]]
				]
			, {orb, Norb}, {l, L-1}, {\[Sigma], f}, {\[Rho], f}];
			V = Table[
				Symbol["V"<>"orb"<>ToString[orb]<>"l"<>ToString[l]<>"\[Sigma]"<>ToString[\[Sigma]]<>"\[Rho]"<>ToString[\[Rho]]]
			, {orb, Norb}, {l, L-1}, {\[Sigma], f}, {\[Rho], f}];
			symbols = {e, V};,
		(* ---------------------------------------------------- *)
			EdMode == "Magnetic",
			e = Table[
					Symbol["e"<>"orb"<>ToString[orb]<>"l"<>ToString[l]<>"\[Sigma]"<>ToString[\[Sigma]]]
			, {orb, Norb}, {l, L-1}, {\[Sigma], f}];
			V = Table[
				Symbol["V"<>"orb"<>ToString[orb]<>"l"<>ToString[l]<>"\[Sigma]"<>ToString[\[Sigma]]]
			, {orb, Norb}, {l, L-1}, {\[Sigma], f}];
			symbols = {e, V};,
		(* -------------------------------------------------- *)
			EdMode == "FullSuperc",
			e = Table[
				ArrayFlatten[{
					{Symbol["e"<>ToString[i]<>ToString[1]] * PauliMatrix[3] + Symbol["\[CapitalDelta]"<>ToString[i]<>ToString[1]] * PauliMatrix[1],
					Symbol["\[CapitalXi]"<>ToString[i]] * PauliMatrix[1]},
					{Symbol["\[CapitalXi]"<>ToString[i]] * PauliMatrix[1], 
					Symbol["e"<>ToString[i]<>ToString[2]] * PauliMatrix[3] + Symbol["\[CapitalDelta]"<>ToString[i]<>ToString[2]] * PauliMatrix[1]}
				}]
			, {i, 1, L-1}];
			V = Table[
				ArrayFlatten[{
					{Symbol["V"<>ToString[i]<>ToString[1]] * PauliMatrix[3],
					0 * PauliMatrix[1]},
					{0 * PauliMatrix[1], 
					Symbol["V"<>ToString[i]<>ToString[2]] * PauliMatrix[3]}
				}]
			, {i, 1, L-1}];
			symbols = {e, V};
		];
		Export[filename, symbols];
	];
	symbols
];

(* implement default symmetries to the symbols representing bath parameters *)
SymmetrizeSymbols[symbols_, EdMode_, OrbitalSymmetry_, ParticleHoleSymmetry_] := Module[
	{e, V, \[CapitalDelta], Norb, Nbath, f},
	If[!OrbitalSymmetry && !ParticleHoleSymmetry, Return[symbols]];
	(* ------------- *)
	Which[
		EdMode == "Normal", 
		e = symbols[[1]]; V = symbols[[2]];
		{Norb, Nbath} = Dimensions[e];
		If[ParticleHoleSymmetry,
			Do[
				If[OddQ[Nbath], 
					e[[orb]][[(Nbath+1)/2]] = 0; 
					Do[
						e[[orb]][[Nbath + 1 - l]] = -e[[orb]][[l]];
						V[[orb]][[Nbath + 1 - l]] = V[[orb]][[l]];
					, {l, (Nbath-1)/2}]
				];
				If[EvenQ[Nbath],
					Do[
						e[[orb]][[Nbath + 1 - l]] = -e[[orb]][[l]];
						V[[orb]][[Nbath + 1 - l]] = V[[orb]][[l]];
					, {l, Nbath/2}]
				];
			, {orb, Norb}]
		];
		If[OrbitalSymmetry, 
			e = ConstantArray[e[[1]], Norb];
			V = ConstantArray[V[[1]], Norb];
		];
		Return[{e,V}];,
	(* ---------------------------------------- *)
		EdMode == "Superc", 
		e = symbols[[1]]; V = symbols[[2]]; \[CapitalDelta] = symbols[[3]];
		{Norb, Nbath} = Dimensions[e];
		If[ParticleHoleSymmetry,
			Do[
				If[OddQ[Nbath], 
					e[[orb]][[(Nbath+1)/2]] = 0; 
					Do[
						e[[orb]][[Nbath + 1 - l]] = -e[[orb]][[l]];
						V[[orb]][[Nbath + 1 - l]] = V[[orb]][[l]];
						\[CapitalDelta][[orb]][[Nbath + 1 - l]] = \[CapitalDelta][[orb]][[l]];
					, {l, (Nbath-1)/2}]
				];
				If[EvenQ[Nbath],
					Do[
						e[[orb]][[Nbath + 1 - l]] = -e[[orb]][[l]];
						V[[orb]][[Nbath + 1 - l]] = V[[orb]][[l]];
						\[CapitalDelta][[orb]][[Nbath + 1 - l]] = \[CapitalDelta][[orb]][[l]];
					, {l, Nbath/2}]
				];
			, {orb, Norb}]
		];
		If[OrbitalSymmetry, 
			e = ConstantArray[e[[1]], Norb];
			V = ConstantArray[V[[1]], Norb];
			\[CapitalDelta] = ConstantArray[\[CapitalDelta][[1]], Norb];
		];
		Return[{e,V,\[CapitalDelta]}];,
	(* ---------------------------------------- *)
		EdMode == "Raman",
		e = symbols[[1]]; V = symbols[[2]];
		{Norb, Nbath, f} = Dimensions[e][[;;3]];
		If[ParticleHoleSymmetry,
			Do[
				If[OddQ[Nbath], 
					Do[
						e[[orb]][[(Nbath+1)/2]][[\[Sigma],\[Sigma]]] = 0; 
						Do[
							If[\[Sigma] != \[Rho], V[[orb]][[(Nbath+1)/2]][[\[Sigma],\[Rho]]] = 0; ]
						, {\[Rho], f}]
					, {\[Sigma], f}];
					Do[
						If[\[Rho] == \[Sigma],
							e[[orb]][[Nbath + 1 - l]][[\[Sigma],\[Sigma]]] = -e[[orb]][[l]][[\[Sigma],\[Sigma]]];
							V[[orb]][[Nbath + 1 - l]][[\[Sigma],\[Sigma]]] = -V[[orb]][[l]][[\[Sigma],\[Sigma]]];,
						(* else *)
							e[[orb]][[Nbath + 1 - l]][[\[Sigma],\[Rho]]] = e[[orb]][[l]][[\[Sigma],\[Rho]]];
							V[[orb]][[Nbath + 1 - l]][[\[Sigma],\[Rho]]] = 0;
							V[[orb]][[l]][[\[Sigma],\[Rho]]] = 0;
						]
					, {l, (Nbath-1)/2}, {\[Sigma], f}, {\[Rho], f}]
				];
				If[EvenQ[Nbath],
					Do[
						If[\[Sigma] == \[Rho],
							e[[orb]][[Nbath + 1 - l]][[\[Sigma],\[Sigma]]] = -e[[orb]][[l]][[\[Sigma],\[Sigma]]];
							V[[orb]][[Nbath + 1 - l]][[\[Sigma],\[Sigma]]] = -V[[orb]][[l]][[\[Sigma],\[Sigma]]];,
						(* else *)
							e[[orb]][[Nbath + 1 - l]][[\[Sigma],\[Rho]]] = e[[orb]][[l]][[\[Sigma],\[Rho]]];
							V[[orb]][[Nbath + 1 - l]][[\[Sigma],\[Rho]]] = 0;
							V[[orb]][[l]][[\[Sigma],\[Rho]]] = 0;
						];
					, {l, Nbath/2}, {\[Sigma], f}, {\[Rho], f}]
				];
			, {orb, Norb}]
		];
		If[OrbitalSymmetry, 
			e = ConstantArray[e[[1]], Norb];
			V = ConstantArray[V[[1]], Norb];
		];
		Return[{e,V}];,
	(* ---------------------------------------- *)
		EdMode == "Magnetic",
		e = symbols[[1]]; V = symbols[[2]];
		{Norb, Nbath, f} = Dimensions[e];
		If[ParticleHoleSymmetry,
			Do[
				If[OddQ[Nbath], 
					Do[
						e[[orb]][[(Nbath+1)/2]][[f+1-\[Sigma]]] = -e[[orb]][[(Nbath+1)/2]][[\[Sigma]]];
						V[[orb]][[(Nbath+1)/2]][[f+1-\[Sigma]]] = V[[orb]][[(Nbath+1)/2]][[\[Sigma]]];		
					, {\[Sigma], f}];
					Do[
						If[\[Sigma] == (f+1)/2, 
							e[[orb]][[Nbath + 1 - l]][[\[Sigma]]] = e[[orb]][[l]][[\[Sigma]]];
							V[[orb]][[Nbath + 1 - l]][[\[Sigma]]] = V[[orb]][[l]][[\[Sigma]]];
							Continue[];
						];
						e[[orb]][[Nbath + 1 - l]][[f+1-\[Sigma]]] = -e[[orb]][[l]][[\[Sigma]]];
						V[[orb]][[Nbath + 1 - l]][[f+1-\[Sigma]]] = V[[orb]][[l]][[\[Sigma]]];
					, {l, (Nbath-1)/2}, {\[Sigma], f}]
				];
				If[EvenQ[Nbath],
					Do[
						If[\[Sigma] == (f+1)/2, 
							e[[orb]][[Nbath + 1 - l]][[\[Sigma]]] = e[[orb]][[l]][[\[Sigma]]];
							V[[orb]][[Nbath + 1 - l]][[\[Sigma]]] = V[[orb]][[l]][[\[Sigma]]];
							Continue[];
						];
						e[[orb]][[Nbath + 1 - l]][[f+1-\[Sigma]]] = -e[[orb]][[l]][[\[Sigma]]];
						V[[orb]][[Nbath + 1 - l]][[f+1-\[Sigma]]] = V[[orb]][[l]][[\[Sigma]]];
					, {l, Nbath/2}, {\[Sigma], f}]
				];
			, {orb, Norb}]
		];
		If[OrbitalSymmetry, 
			e = ConstantArray[e[[1]], Norb];
			V = ConstantArray[V[[1]], Norb];
		];
		Return[{e,V}];
	]
];

(* generate flat list with the independent symbols *)
IndependentSymbols[symbols_] := DeleteDuplicates[
	DeleteDuplicates[
		DeleteCases[Flatten[symbols], 0]
	], #1 == -#2&
];

(* list of indexes of the position of independent symbols in the tensor symbol *)
IndependentSymbolsIndexes[symbols_, independentsymbols_] := First[#] &/@ (Position[symbols, #] &/@ independentsymbols);

(* extract independent bath parameters depending on EdMode *)
(*TakeIndependentParameters[L_, f_, Norb_, \[Sigma]_, orb_, BathParameters_, EdMode_] :=
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
	];*)
(* new version *)
TakeIndependentParameters[BathParameters_, independentsymbolindexes_] := MapApply[
	Part[BathParameters, ##]&, independentsymbolindexes
];

(* correctly reshape the flat list of independent bath parameters *)
(*ReshapeBathParameters[L_, f_, Norb_, IndependentParameters_, OrbitalSymmetry_, EdMode_] := 
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
		{(* e *)
			ConstantArray[(
					SymmetricMatrixFromArray[#,f] &/@ Partition[IndependentParameters, f*(f+1)/2]
				)[[ ;; L-1]]
			, Norb],
		 (* V *)
		 	ConstantArray[(
					SymmetricMatrixFromArray[#,f] &/@ Partition[IndependentParameters, f*(f+1)/2]
				)[[-(L-1) ;; ]]
			, Norb]
		},
	(* ----------------------------------------- *)
		EdMode == "Raman" && !OrbitalSymmetry,
		{(* e *)
			Table[(
					SymmetricMatrixFromArray[#,f] &/@ Partition[IndependentParameters[[orb]], f*(f+1)/2]
				)[[ ;; L-1]]
			, {orb, Norb}],
		 (* V *)
		 	Table[(
					SymmetricMatrixFromArray[#,f] &/@ Partition[IndependentParameters[[orb]], f*(f+1)/2]
				)[[-(L-1) ;; ]]
			, {orb, Norb}]
		},	
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
	];*)
(* new version *)
ReshapeBathParameters[symbols_, independentsymbols_, IndependentBathParameters_] := Module[
	{rules = Thread[independentsymbols -> IndependentBathParameters]},
	symbols /. rules
];

End[]

EndPackage[]
