(* ::Package:: *)

BeginPackage["Hamiltonian`", {"Sectors`"}]


HNonlocal::usage = "HNonlocal[L, f, Norb, Sectors, EdMode]"

HNonlocalRaman::usage = "HNonlocalRaman[L, f, Norb, Sectors, EdMode]"

HLocal::usage = "HLocal[L, f, Norb, Sectors, EdMode] "

GetHamiltonian::usage = "GetHamiltonian[L_, f_, Norb_, Sectors_, LoadHamiltonianQ_, HnonlocFile_, HlocFile_, EdMode_]"

HImp::usage = "HImp[Norb_, HnonlocBlocks_, HlocBlocks_, BathParameters_, InteractionParameters_, EdMode_]"

Hnonint::usage = "Hnonint[L, f, Norb, Sectors, EdMode]. Optional arguments: RealPBC -> True (default)/ False; SyntheticPBC -> True / False (default), RealPhase -> {0}, SyntheticPhase -> 0  "


Begin["Private`"]

Print["Package Hamiltonian` loaded successfully."];

(* gives right neighbor of a site j in a closed (open) chain of L elements *)
Neighbor = Compile[{
	{L,_Integer}, {j,_Integer}
	},
	Mod[j,L]+1,
	CompilationTarget->"C"
];


(*          BUILD THE IMPURITY HAMILTONIAN           *)
(* Non-local Hamiltonian blocks *)
HNonlocal[L_, f_, Norb_, Sectors_, EdMode_, OptionsPattern[]] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H = {};
	Do[
		Hsector = {};
		dim = Length[\[Psi]];
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Hblock = SparseArray[{}, {dim,dim}];
			Which[
				flag == "Bath" && EdMode != "Raman",
				num = n[L, f, Norb, j, \[Sigma], orb, \[Psi]];(*local density*)
				Hblock += SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Hopping" && EdMode != "Raman",
				\[Psi]1 = HopSelect[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]];
				If[Length[\[Psi]1] != 0, 
					\[Chi] = Hop[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]1];
					rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
					\[CapitalSigma] = CCSign[L, f, {1,j}, {\[Sigma],\[Sigma]}, {orb,orb}, \[Psi]1];
					Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					Hblock = Hblock + Hblock\[ConjugateTranspose];
				];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Superc" && (EdMode == "Superc" || EdMode == "FullSuperc"),
				If[\[Sigma] == 1,
					\[Psi]1 = CreatePairSelect[L, f, j, j, 1, 2, orb, orb, \[Psi]];
					If[Length[\[Psi]1] != 0, 
						\[Chi] = CreatePair[L, f, j, j, 1, 2, orb, orb, \[Psi]1];
						rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						\[CapitalSigma] = CCSign[L, f, {j,j}, {1,2}, {orb,orb}, \[Psi]1];
						Hblock += SparseArray[pos->\[CapitalSigma], {dim,dim}];
						Hblock = Hblock + Hblock\[ConjugateTranspose];
					];
					AppendTo[Hsector, Hblock];
				];,
			(* --------------------------------------------------------------- *)
				flag == "InterorbSuperc" && (EdMode == "InterorbSuperc" || EdMode == "FullSuperc"),
				If[\[Sigma] == 1 && orb == 1, (* <- do NOT provide a block for every spin and orbital index! *)
				Do[
					If[
						orb2 != orb1,
						\[Psi]1 = CreatePairSelect[L, f, j, j, 1, 2, orb1, orb2, \[Psi]];
						If[Length[\[Psi]1] != 0,
							\[Chi] = CreatePair[L, f, j, j, 1, 2, orb1, orb2, \[Psi]1];
							rows = \[Chi]/.dispatch;(* *)cols = \[Psi]1/.dispatch;(* *)pos = {rows,cols}\[Transpose];
							\[CapitalSigma] = CCSign[L, f, {j,j}, {1,2}, {orb1,orb2}, \[Psi]1];
							If[orb2 < orb1, \[CapitalSigma] = -\[CapitalSigma]]; (* ?? *)
							Hblock += SparseArray[pos -> \[CapitalSigma], {dim,dim}];
						];
					]
				, {orb1, Norb}, {orb2, Norb}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];
				];
			];
		,{flag, {"Bath","Hopping","Superc","InterorbSuperc"}}, {orb,1,Norb}, {\[Sigma],1,f}, {j,OptionValue["Nimp"]+1,L}];
		AppendTo[H, Hsector];
	,{\[Psi], Sectors}];
	H
];
Options[HNonlocal] = {"Nimp" -> 1};

(* non local part of impurity Hamiltonian for Raman processes *)
HNonlocalRaman[L_, f_, Norb_, Sectors_, EdMode_, OptionsPattern[]] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H = {};
	Do[
		Hsector = {};
		dim = Length[\[Psi]];
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Hblock = SparseArray[{}, {dim,dim}];
			Which[
				flag == "Bath" && EdMode == "Raman",
				(*Print["flag=",flag,". orb=", orb,". {\[Rho],\[Sigma]}=",{\[Rho],\[Sigma]},". j=",j];*)
				\[Psi]1 = HopSelect[L, f, j, j, \[Rho], \[Sigma], orb, orb, \[Psi]];
				If[Length[\[Psi]1] != 0, 
					\[Chi] = Hop[L, f, j, j, \[Rho], \[Sigma], orb, orb, \[Psi]1];
					rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
					\[CapitalSigma] = CCSign[L, f, {j,j}, {\[Rho],\[Sigma]}, {orb,orb}, \[Psi]1];
					Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					Hblock = Hblock + Hblock\[ConjugateTranspose];
				];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Hopping" && EdMode == "Raman",
				(*Print["flag=",flag,". orb=", orb,". {\[Rho],\[Sigma]}=",{\[Rho],\[Sigma]},". j=",j];*)
				\[Psi]1 = HopSelect[L, f, 1, j, \[Rho], \[Sigma], orb, orb, \[Psi]];
				If[Length[\[Psi]1] != 0, 
					\[Chi] = Hop[L, f, 1, j, \[Rho], \[Sigma], orb, orb, \[Psi]1];
					rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
					\[CapitalSigma] = CCSign[L, f, {1,j}, {\[Rho],\[Sigma]}, {orb,orb}, \[Psi]1];
					Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					Hblock = Hblock + Hblock\[ConjugateTranspose];
				];
				AppendTo[Hsector, Hblock];
			];
		,{flag, {"Bath", "Hopping"}}, {orb, 1, Norb}, {j, OptionValue["Nimp"]+1, L}, {\[Rho], 1, f}, {\[Sigma], 1, f}];
		AppendTo[H, Hsector];
	,{\[Psi], Sectors}];
	H
];
Options[HNonlocalRaman] = {"Nimp" -> 1};

(* Local Hamiltonian blocks *)
HLocal[L_, f_, Norb_, Sectors_, EdMode_, OptionsPattern[]] := Module[
	{H,Hsector,Hblock,rules,dispatch,num,dim,\[Psi]1,\[Chi],rows,cols,pos,\[CapitalSigma]},
	H = {};
	Do[
		Hsector = {};
		dim = Length[\[Psi]];
		Hsector = {};
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Which[
				flag == "Hubbard",
				Do[
					Hblock = SparseArray[{},{dim,dim}];
					Do[
						num = Sum[
							n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
						, {\[Sigma], 1, f}];
						Hblock += SparseArray @ DiagonalMatrix[0.5*num*(num-1)];
					, {j, 1, OptionValue["Nimp"]}];
					AppendTo[Hsector, Hblock];
				, {orb, 1, Norb}],
			(* ---------------------------------- *)
				flag == "Interorb_Hubbard_Opposite_Spin" && Norb > 1,
				num = Sum[
					If[orbA != orbB,
						n[L,f,Norb,j,1,orbA,\[Psi]]*n[L,f,Norb,j,2,orbB,\[Psi]],
					(*else*)
						0]
				,{orbA, Norb}, {orbB, Norb}, {j, 1, OptionValue["Nimp"]}];
				Hblock = SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];,
			(* ---------------------------------- *)
				flag == "Interorb_Hubbard_Same_Spin" && Norb > 1,
				num = Sum[
					n[L, f, Norb, j, \[Sigma], orb, \[Psi]] * n[L, f, Norb, j, \[Sigma], orb+1, \[Psi]]
				,{\[Sigma], f}, {orb, Norb-1}, {j, 1, OptionValue["Nimp"]}];
				Hblock = SparseArray @ DiagonalMatrix[num];
				AppendTo[Hsector,Hblock];,
			(* ---------------------------------- *)
				(*flag == "Pair_Hopping" && Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "Superc" || EdMode == "FullSuperc"),
				Hblock = SparseArray[{},{dim,dim}];
				Do[
					If[orb2 > orb1,
						\[Psi]1 = PairHoppingSelect[L, f, 1, orb1, orb2, \[Psi]];
						If[Length[\[Psi]1] == 0, Continue[];];
						\[Chi] = PairHopping[L, f, 1, orb1, orb2, \[Psi]1];
						rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						\[CapitalSigma] = CCSign[L, f, {j,j,j,j}, {1,2,2,1}, {orb1,orb1,orb2,orb2}, \[Psi]1];
						Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					];
				,{orb1,1,Norb}, {orb2,1,Norb}, {j,1,OptionValue["Nimp"]}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector,Hblock];*)
				flag == "Pair_Hopping" && Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "Superc" || EdMode == "FullSuperc"),
				Hblock = SparseArray[{}, {dim,dim}];
				Do[
					If[orb2 >= orb1,
						\[Psi]1 = PairHoppingSelect[L, f, 1, orb1, orb2, \[Psi]];
						If[Length[\[Psi]1] == 0, AppendTo[Hsector, Hblock]; Continue[];];
						\[Chi] = PairHopping[L, f, 1, orb1, orb2, \[Psi]1];
						rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						(* adg_up adg_dw b_dw b_up = - adg_up adg_dw b_up b_dw <---- with this choice the operators do not leap over themselves when applied to a state *)
						\[CapitalSigma] = - CCSign[L, f, {j,j,j,j}, {1,2,1,2}, {orb1,orb1,orb2,orb2}, \[Psi]1];
						Hblock = SparseArray[pos -> \[CapitalSigma], {dim,dim}];
						If[orb2 > orb1, Hblock = Hblock + Hblock\[ConjugateTranspose]; ];
						AppendTo[Hsector, Hblock];
					];
				,{orb1,1,Norb}, {orb2,1,Norb}, {j,1,OptionValue["Nimp"]}];,
			(* ----------------------------------- *)
				flag == "Spin_Exchange" && Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "Raman" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc"),
				Hblock = SparseArray[{},{dim,dim}];
				Do[
					If[orb2 > orb1,
						\[Psi]1 = SpinExchangeSelect[L, f, 1, orb1, orb2, \[Psi]];
						If[Length[\[Psi]1] == 0, Continue[];];
						\[Chi] = SpinExchange[L, f, 1, orb1, orb2, \[Psi]1];
						rows = \[Chi]/.dispatch;(* *)cols = \[Psi]1/.dispatch;(* *)pos = {rows,cols}\[Transpose];
						\[CapitalSigma] = CCSign[L, f, {j,j,j,j}, {1,2,1,2}, {orb1,orb1,orb2,orb2}, \[Psi]1];
						Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					];
				,{orb1,1,Norb}, {orb2,1,Norb}, {j,1,OptionValue["Nimp"]}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];,
			(* ---------------------------------- *)
				flag == "Energy_Shift",
				num = Sum[
					n[L,f,Norb,j,\[Sigma],orb,\[Psi]]
				,{\[Sigma],1,f}, {orb,1,Norb}, {j,1,OptionValue["Nimp"]}];
				Hblock = SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];,
			(* ----------------------------------- *)
				flag == "Magnetic_Field",
				Do[
					num = Sum[
						n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
					, {j, 1, OptionValue["Nimp"]}];
					Hblock = SparseArray @ DiagonalMatrix[num];
					AppendTo[Hsector, Hblock];
				, {orb, Norb}, {\[Sigma], f}],
			(* ----------------------------------- *)
				flag == "Crystal_Field", (* different on site energy to different orbitals *) 
				Do[
					num = Sum[
						n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
					, {j, 1, OptionValue[Nimp]}, {\[Sigma], f}];
					Hblock = SparseArray @ DiagonalMatrix[num];
					AppendTo[Hsector, Hblock];
				, {orb, Norb}]
			];
		,{flag, {"Magnetic_Field","Crystal_Field","Hubbard","Interorb_Hubbard_Opposite_Spin","Interorb_Hubbard_Same_Spin","Pair_Hopping","Spin_Exchange","Energy_Shift"}}];
		AppendTo[H, Hsector];
	,{\[Psi],Sectors}];
	H
];
Options[HLocal] = {"Nimp" -> 1};

(* Get the Hamiltonian structure once for all *)
GetHamiltonian[L_, f_, Norb_, Nimp_, Sectors_, LoadHamiltonianQ_, HnonlocFile_, HlocFile_, EdMode_] := Module[
	{HnonlocBlocks, HlocBlocks},
	If[LoadHamiltonianQ,
		Print["Getting Hamiltonians from file"];
		HnonlocBlocks = Import[HnonlocFile];
		HlocBlocks = Import[HlocFile];
		Print["Done! Let's get started! \n"],
	(* else *)
		Print["Computing Hamiltonians..."];
		Print["Time: ", First @ AbsoluteTiming[
			If[EdMode != "Raman",
				HnonlocBlocks = HNonlocal[L, f, Norb, Sectors, EdMode, "Nimp" -> Nimp];,
			(* else, EdMode == "Raman" *)
				HnonlocBlocks = HNonlocalRaman[L, f, Norb, Sectors, EdMode, "Nimp" -> Nimp];
			];
			Export[HnonlocFile, HnonlocBlocks];
			HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode, "Nimp" -> Nimp];
			Export[HlocFile, HlocBlocks];
		]," sec."];
		Print["Done! Let's get started! \n"];
	];
	{HnonlocBlocks, HlocBlocks}
];

(* build the impurity Hamiltonian from the local and nonlocal blocks and the respective parameters *)
HImp[Norb_, HnonlocBlocks_, HlocBlocks_, BathParameters_, InteractionParameters_, EdMode_] := Module[
	{EffectiveInteractionParameters, FlatBathParameters = Flatten[BathParameters], Hloc, Hnonloc},
	(* general structure of interaction parameters: {\[Delta], U, Ust, Usec, Jph, Jse, - \[Mu] + shift} *)
	EffectiveInteractionParameters = Which[
		Norb == 1,
		(* with 1 orbital, delete Jse, Jph, Usec, Ust, at positions -2, -3, -4, -5 *)
		Flatten[
			Delete[InteractionParameters, {{-2},{-3},{-4},{-5}}]
		],
	(* --------------------------------------------------- *)
		Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc"),
		(* this is the most complicated scenario, you don't remove anything here *)
		Flatten[InteractionParameters],
	(* --------------------------------------------------- *)
		Norb > 1 && EdMode == "Superc", 
		(* pair hopping is possible, but spin exchange is not *)
		Flatten[
			Delete[InteractionParameters, -2]
		],
	(* --------------------------------------------------- *)
		Norb > 1 && EdMode == "Raman",
		(* pair hopping is not possible (it changes the number of particles per orbital), but spin exchange is ok *)
		Flatten[
			Delete[InteractionParameters, -3]
		],
	(* --------------------------------------------------- *)
		Norb > 1 && EdMode == "Normal",
		(* pair hopping and spin exchange are not possible *)
		Flatten[
			Delete[InteractionParameters, {{-2},{-3}}]
		]
	];
	(* build non local Hamiltonian -> avoid Dot[] as it returns dense array *)
	Hnonloc = SparseArray[#]&/@(
		Sum[
			FlatBathParameters[[i]]*#[[i]],
		{i, 1, Length@FlatBathParameters}] &/@ HnonlocBlocks
	);
	(* build local Hamiltonian -> avoid Dot[] as it returns dense array *)
	Hloc = SparseArray[#]&/@(
		Sum[
			EffectiveInteractionParameters[[i]]*#[[i]],
		{i, 1, Length@EffectiveInteractionParameters}] &/@ HlocBlocks
	);
	Hnonloc + Hloc
];



(*       BUILD CHAIN HAMILTONIAN        *)
Hnonint[L_, f_, Norb_, Sectors_, EdMode_, OptionsPattern[]] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H = {};
	Do[
		Hsector = {};
		dim = Length[\[Psi]];
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Hblock = SparseArray[{}, {dim,dim}];
			Which[
				flag == "Potential",
				num = n[L, f, Norb, j, \[Sigma], orb, \[Psi]];(*local density*)
				Hblock += SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Hopping",
				If[j == L && !OptionValue[RealPBC], Continue[];];
				\[Psi]1 = HopSelect[L, f, j, Neighbor[L,j], \[Sigma], \[Sigma], orb, orb, \[Psi]];
				If[Length[\[Psi]1] == 0, AppendTo[Hsector, Hblock]; Continue[];];
				\[Chi] = Hop[L, f, j, Neighbor[L,j], \[Sigma], \[Sigma], orb, orb, \[Psi]1];
				rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
				\[CapitalSigma] = CCSign[L,f,{j,Neighbor[L,j]},{\[Sigma],\[Sigma]},{orb,orb},\[Psi]1];
				If[j == L, \[CapitalSigma] = -\[CapitalSigma]];(* if j=L, we are applying cdg_L c_1. When moving cdg_L before position L, it jumps over c_1 and the sign changes! *)
				If[
					OptionValue[RealPhase] == {0},
					Hblock += SparseArray[pos -> \[CapitalSigma], {dim,dim}];,
				(* else *)
					Hblock += SparseArray[pos -> \[CapitalSigma]*Exp[I*OptionValue[RealPhase][[\[Sigma]]]], {dim,dim}]
				];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Raman" && EdMode == "Raman",
				If[\[Sigma] == f && !OptionValue[SyntheticPBC], Continue[];];
				\[Psi]1 = HopSelect[L, f, j, j, \[Sigma], Neighbor[f,\[Sigma]], orb, orb, \[Psi]];
				If[Length[\[Psi]1] == 0, AppendTo[Hsector, Hblock]; Continue[];];
				\[Chi] = Hop[L, f, j, j, \[Sigma], Neighbor[f,\[Sigma]], orb, orb, \[Psi]1];
				rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
				\[CapitalSigma] = CCSign[L,f,{j,j},{\[Sigma],Neighbor[f,\[Sigma]]},{orb,orb},\[Psi]1];
				If[\[Sigma] == f, \[CapitalSigma] = -\[CapitalSigma]];(* if \[Sigma]=f, we are applying cdg_f c_1. When moving cdg_f before position f, it jumps over c_1 and the sign changes! *)
				Hblock += SparseArray[pos -> \[CapitalSigma]*Exp[I*OptionValue[SyntheticPhase]*j], {dim,dim}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];	
			];
		,{flag, {"Potential","Hopping","Raman"}}, {orb,1,Norb}, {\[Sigma],1,f}, {j,1,L}];
		AppendTo[H, Hsector];
	,{\[Psi], Sectors}];
	H
];
Options[Hnonint] = {RealPBC -> True, SyntheticPBC -> False, RealPhase -> {0}, SyntheticPhase -> 0};

End[]

EndPackage[]
