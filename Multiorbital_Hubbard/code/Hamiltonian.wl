(* ::Package:: *)

BeginPackage["Hamiltonian`", {"Sectors`"}]


HNonlocalNormal::usage = "HNonlocalNormal[L_, f_, Norb_, Sectors_]"

HLocalNormal::usage = "HLocalNormal[L_, f_, Norb_, Sectors_, EdMode_, InteractionParameters_]"

HLocal::usage = "HLocal[L_, f_, Norb_, Sectors_, EdMode_, InteractionParameters_, OptionsPattern[]]"

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
(* tensor of dimension Nsectors x 2 x Norb x Nbath whose objects are SparseArrays *)
HNonlocalNormal[L_, f_, Norb_, Sectors_] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = Last @ Last @ Reap[
				Do[
					Which[
						flag == "Bath",
						num = Sum[n[L, f, Norb, j, \[Sigma], orb, \[Psi]] , {\[Sigma], f}];(*local density*)
						Hblock = SparseArray @ DiagonalMatrix[num];
						Sow[Hblock, "Hblock"];,
					(* --------------------------------------------------------------- *)
						flag == "Hopping",
						Hblock = SparseArray[{}, {dim,dim}];
						Do[
							\[Psi]1 = HopSelect[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]];
							If[Length[\[Psi]1] == 0, Continue[]; ]; 
							\[Chi] = Hop[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]1];
							rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
							\[CapitalSigma] = CCSign[L, f, {1,j}, {\[Sigma],\[Sigma]}, {orb,orb}, \[Psi]1];
							Hblock += SparseArray[pos -> \[CapitalSigma], {dim, dim}];
						, {\[Sigma], f}];					
						Sow[Hblock + Hblock\[ConjugateTranspose], "Hblock"];
					];
				, {flag, {"Bath", "Hopping"}}, {orb, Norb}, {j, 2, L}]
			, "Hblock"];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}]
	, "Hsector"]
];
(* tensor of dimension Nsectors x 3 x Norb x Nbath whose objects are SparseArrays *)
HNonlocalSuperc[L_, f_, Norb_, Sectors_] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = Last @ Last @ Reap[
				Do[
					Which[
						flag == "Bath",
						num = Sum[n[L, f, Norb, j, \[Sigma], orb, \[Psi]] , {\[Sigma], f}];(*local density*)
						Hblock = SparseArray @ DiagonalMatrix[num];
						Sow[Hblock, "Hblock"];,
					(* --------------------------------------------------------------- *)
						flag == "Hopping",
						Hblock = SparseArray[{}, {dim,dim}];
						Do[
							\[Psi]1 = HopSelect[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]];
							If[Length[\[Psi]1] == 0, Continue[]; ]; 
							\[Chi] = Hop[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]1];
							rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
							\[CapitalSigma] = CCSign[L, f, {1,j}, {\[Sigma],\[Sigma]}, {orb,orb}, \[Psi]1];
							Hblock += SparseArray[pos -> \[CapitalSigma], {dim, dim}];
						, {\[Sigma], f}];					
						Sow[Hblock + Hblock\[ConjugateTranspose], "Hblock"];,
					(* --------------------------------------------------------------- *)
						flag == "Superc",
						Hblock = SparseArray[{}, {dim,dim}];
						\[Psi]1 = CreatePairSelect[L, f, j, j, 1, 2, orb, orb, \[Psi]];
						If[Length[\[Psi]1] != 0, 
							\[Chi] = CreatePair[L, f, j, j, 1, 2, orb, orb, \[Psi]1];
							rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
							\[CapitalSigma] = CCSign[L, f, {j,j}, {1,2}, {orb,orb}, \[Psi]1];
							Hblock += SparseArray[pos -> \[CapitalSigma], {dim, dim}];
						];
						Sow[Hblock + Hblock\[ConjugateTranspose], "Hblock"];
					];
				, {flag, {"Bath", "Hopping", "Superc"}}, {orb, 1, Norb}, {j, 2, L}]
			, "Hblock"];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}]
	, "Hsector"]
];
(* non local part of impurity Hamiltonian for Raman processes *)
HNonlocalRaman[L_, f_, Norb_, Sectors_] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = Last @ Last @ Reap[
				Do[		
					Hblock = SparseArray[{}, {dim,dim}];			
					Which[
						flag == "Bath",
						If[\[Rho] == \[Sigma],
							num = n[L, f, Norb, j, \[Sigma], orb, \[Psi]];(* local density*)
							Hblock += SparseArray[Band[{1,1}] -> num];,
						(* else \[Rho] != \[Sigma] *)
							(*If[OptionValue["OnlyDiagonal"], Continue[]; ]; 
							If[OptionValue["OnlyTriangular"] && \[Sigma] < \[Rho], Continue[]; ];*)
							\[Psi]1 = HopSelect[L, f, j, j, \[Rho], \[Sigma], orb, orb, \[Psi]];
							If[Length[\[Psi]1] != 0, 
								\[Chi] = Hop[L, f, j, j, \[Rho], \[Sigma], orb, orb, \[Psi]1];
								rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
								\[CapitalSigma] = CCSign[L, f, {j,j}, {\[Rho],\[Sigma]}, {orb,orb}, \[Psi]1];
								If[Index[L, f, Norb, j, \[Rho], orb] > Index[L, f, Norb, j, \[Sigma], orb], 
									\[CapitalSigma] = -\[CapitalSigma]
								];(* if \[Rho] > \[Sigma], we are applying e.g. cdg_2 c_1. When moving cdg_2 before position 2, it jumps over c_1 and the sign changes! *)
								Hblock += SparseArray[pos -> \[CapitalSigma], {dim, dim}];
								(*If[OptionValue["OnlyTriangular"], Hblock = Hblock + Hblock\[ConjugateTranspose]; ];*)
							];
						];
						Sow[Hblock, "Hblock"];,
					(* --------------------------------------------------------------- *)
						flag == "Hopping",
						(*Print["flag=",flag,". orb=", orb,". {\[Rho],\[Sigma]}=",{\[Rho],\[Sigma]},". j=",j];*)
						(*If[OptionValue["OnlyDiagonal"] && \[Sigma] != \[Rho], Continue[]; ]; *)
						\[Psi]1 = HopSelect[L, f, 1, j, \[Rho], \[Sigma], orb, orb, \[Psi]];
						If[Length[\[Psi]1] != 0, 
							\[Chi] = Hop[L, f, 1, j, \[Rho], \[Sigma], orb, orb, \[Psi]1];
							rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
							\[CapitalSigma] = CCSign[L, f, {1,j}, {\[Rho],\[Sigma]}, {orb,orb}, \[Psi]1];
							If[Index[L, f, Norb, 1, \[Rho], orb] > Index[L, f, Norb, j, \[Sigma], orb], 
								\[CapitalSigma] = -\[CapitalSigma]
							];(* if \[Rho] > \[Sigma], we are applying e.g. cdg_2 c_1. When moving cdg_2 before position 2, it jumps over c_1 and the sign changes! *)
							Hblock += SparseArray[pos -> \[CapitalSigma], {dim, dim}];
							Hblock += Hblock\[ConjugateTranspose];
						];
						Sow[Hblock, "Hblock"];
					];
				, {flag, {"Bath", "Hopping"}}, {orb, Norb}, {j, 2, L}, {\[Rho], f}, {\[Sigma], f}]
			, "Hblock"];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"]
];
(* non local part of impurity Hamiltonian for Raman processes *)
HNonlocalMagnetic[L_, f_, Norb_, Sectors_] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = Last @ Last @ Reap[
				Do[		
					Hblock = SparseArray[{}, {dim,dim}];			
					Which[
						flag == "Bath",
						num = n[L, f, Norb, j, \[Sigma], orb, \[Psi]];(* local density*)
						Hblock += SparseArray[Band[{1,1}] -> num];
						Sow[Hblock, "Hblock"];,
					(* --------------------------------------------------------------- *)
						flag == "Hopping",
						(*Print["flag=",flag,". orb=", orb,". {\[Rho],\[Sigma]}=",{\[Rho],\[Sigma]},". j=",j];*)
						(*If[OptionValue["OnlyDiagonal"] && \[Sigma] != \[Rho], Continue[]; ]; *)
						\[Psi]1 = HopSelect[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]];
						If[Length[\[Psi]1] != 0, 
							\[Chi] = Hop[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]1];
							rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
							\[CapitalSigma] = CCSign[L, f, {1,j}, {\[Sigma],\[Sigma]}, {orb,orb}, \[Psi]1];
							If[Index[L, f, Norb, 1, \[Sigma], orb] > Index[L, f, Norb, j, \[Sigma], orb], 
								\[CapitalSigma] = -\[CapitalSigma]
							];(* if \[Rho] > \[Sigma], we are applying e.g. cdg_2 c_1. When moving cdg_2 before position 2, it jumps over c_1 and the sign changes! *)
							Hblock += SparseArray[pos -> \[CapitalSigma], {dim, dim}];
							Hblock += Hblock\[ConjugateTranspose];
						];
						Sow[Hblock, "Hblock"];
					];
				, {flag, {"Bath", "Hopping"}}, {orb, Norb}, {j, 2, L}, {\[Sigma], f}]
			, "Hblock"];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"]
];

(* Local Hamiltonian blocks *)
(*HLocal[L_, f_, Norb_, Sectors_, EdMode_, OptionsPattern[]] := Module[
	{H,Hsector,Hblock,rules,dispatch,num,dim,\[Psi]1,\[Chi],rows,cols,pos,\[CapitalSigma]},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = Last @ Last @ Reap[
				Do[
					Which[
						flag == "Hubbard",
						Do[
							Hblock = SparseArray[{},{dim,dim}];
							Do[
								num = Sum[
									n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
								, {\[Sigma], f}];
								Hblock += SparseArray @ DiagonalMatrix[0.5*num*(num-1)];
							, {j, 1, OptionValue["Nimp"]}];
							Sow[Hblock, "Hblock"];
						, {orb, Norb}],
					(* ---------------------------------- *)
						flag == "Interorb_Hubbard_Opposite_Spin" && Norb > 1,
						num = Sum[
							If[orbA != orbB,
								n[L,f,Norb,j,1,orbA,\[Psi]]*n[L,f,Norb,j,2,orbB,\[Psi]],
							(*else*)
								0]
						,{orbA, Norb}, {orbB, Norb}, {j, 1, OptionValue["Nimp"]}];
						Hblock = SparseArray @ DiagonalMatrix[num];
						Sow[Hblock, "Hblock"];,
					(* ---------------------------------- *)
						flag == "Interorb_Hubbard_Same_Spin" && Norb > 1,
						num = Sum[
							n[L, f, Norb, j, \[Sigma], orb, \[Psi]] * n[L, f, Norb, j, \[Sigma], orb+1, \[Psi]]
						,{\[Sigma], f}, {orb, Norb-1}, {j, 1, OptionValue["Nimp"]}];
						Hblock = SparseArray @ DiagonalMatrix[num];
						Sow[Hblock, "Hblock"];,
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
								Sow[Hblock, "Hblock"];
							];
						,{orb1, Norb}, {orb2, Norb}, {j,1,OptionValue["Nimp"]}];,
					(* ----------------------------------- *)
						flag == "Spin_Exchange" && Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "Raman" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc"),
						Hblock = SparseArray[{},{dim,dim}];
						Do[
							\[Psi]1 = SpinExchangeSelect[L, f, 1, orb1, orb2, \[Psi]];
							If[Length[\[Psi]1] == 0, Continue[];];
							\[Chi] = SpinExchange[L, f, 1, orb1, orb2, \[Psi]1];
							rows = \[Chi]/.dispatch;(* *)cols = \[Psi]1/.dispatch;(* *)pos = {rows,cols}\[Transpose];
							\[CapitalSigma] = -CCSign[L, f, {j,j,j,j}, {1,2,2,1}, {orb1,orb1,orb2,orb2}, \[Psi]1];
							Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
						,{orb1, 1, Norb}, {orb2, orb1+1, Norb}, {j,1,OptionValue["Nimp"]}];
						Hblock = Hblock + Hblock\[ConjugateTranspose];
						Sow[Hblock, "Hblock"];,
					(* ---------------------------------- *)
						flag == "Energy_Shift",
						num = Sum[
							n[L,f,Norb,j,\[Sigma],orb,\[Psi]]
						, {\[Sigma], f}, {orb, Norb}, {j, 1, OptionValue["Nimp"]}];
						Hblock = SparseArray @ DiagonalMatrix[num];
						Sow[Hblock, "Hblock"];,
					(* ----------------------------------- *)
						flag == "Magnetic_Field", (* this should be deprecated in favor of the more general Raman field *) 
						Do[
							num = Sum[
								n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
							, {j, 1, OptionValue["Nimp"]}];
							Hblock = SparseArray @ DiagonalMatrix[num];
							Sow[Hblock, "Hblock"];
						, {orb, Norb}, {\[Sigma], f}],
					(* ----------------------------------- *)
						flag == "Raman_Field", (* this is coupled to a flattened list of Norb uppertriangularized Raman fxf matrices *)
						Do[
							If[\[Rho] == \[Sigma], (* diagonal elements of Raman matrix represent a magnetic field *)
								num = Sum[
									n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
								, {j, 1, OptionValue["Nimp"]}];
								Hblock = SparseArray @ DiagonalMatrix[num];,
							(* else, if \[Rho] > \[Sigma], off diagonal elements of Raman matrix *)
								Hblock = SparseArray[{},{dim,dim}];
								\[Psi]1 = HopSelect[L, f, 1, 1, \[Rho], \[Sigma], orb, orb, \[Psi]];
								If[Length[\[Psi]1] != 0, 
									\[Chi] = Hop[L, f, 1, 1, \[Rho], \[Sigma], orb, orb, \[Psi]1];
									rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
									\[CapitalSigma] = CCSign[L, f, {1,1}, {\[Rho],\[Sigma]}, {orb,orb}, \[Psi]1];
									If[Index[L, f, Norb, 1, \[Rho], orb] > Index[L, f, Norb, 1, \[Sigma], orb], 
										\[CapitalSigma] = -\[CapitalSigma]
									];(* if \[Rho] > \[Sigma], we are applying e.g. cdg_2 c_1. When moving cdg_2 before position 2, it jumps over c_1 and the sign changes! *)
									Hblock = SparseArray[pos->\[CapitalSigma],{dim,dim}];
									Hblock = Hblock + Hblock\[ConjugateTranspose];
								];
							];
							Sow[Hblock, "Hblock"];
						, {orb, Norb}, {\[Rho], f}, {\[Sigma], \[Rho], f}],
					(* ----------------------------------- *)
						flag == "Crystal_Field", (* different on site energy to different orbitals *) 
						Do[
							num = Sum[
								n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
							, {j, 1, OptionValue[Nimp]}, {\[Sigma], f}];
							Hblock = SparseArray @ DiagonalMatrix[num];
							Sow[Hblock, "Hblock"];
						, {orb, Norb}]
					];
				,{flag, {"Raman_Field","Crystal_Field","Hubbard","Interorb_Hubbard_Opposite_Spin","Interorb_Hubbard_Same_Spin","Pair_Hopping","Spin_Exchange","Energy_Shift"}}];
			, "Hblock"];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"]
];
Options[HLocal] = {"Nimp" -> 1};*)

(* only 1. Hubbard interaction; 2. Ust; 3. Usec; 4. chemical potential; 5. Crystal field *)
HLocalNormal[L_, f_, Norb_, Sectors_, EdMode_, InteractionParameters_, OptionsPattern[]] := Module[{
	U = InteractionParameters[[1]],
	Ust = InteractionParameters[[2]],
	Usec = InteractionParameters[[3]],
	\[Mu]eff = InteractionParameters[[4]],
	\[Delta] = InteractionParameters[[5]],
	H,Hsector,Hblock,rules,dispatch,num,norb,dim,\[Psi]1,\[Chi],rows,cols,pos,\[CapitalSigma]
	},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = SparseArray[{}, {dim, dim}];
			(* exception: if all parameters vanish, avoid calculations *)
			If[# == ConstantArray[0, Length[#]],
				Sow[Hsector, "Hsector"];
				Continue[];,
			(* else, compute n_j,\[Sigma],orb on all the states of the sector *)
				num = Table[
					n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
				, {j, 1, OptionValue["Nimp"]}, {\[Sigma], f}, {orb, Norb}];
			] &@ Flatten[InteractionParameters];
			(* orbital-wise Hubbard term *)
			Do[
				If[U[[orb]] == 0, Continue[]; ];
				norb = Sum[num[[j,\[Sigma],orb]], {\[Sigma],f}];
				Hsector += 0.5*U[[orb]]* SparseArray @ DiagonalMatrix[
					norb * (norb - 1)
				];		
			, {orb, Norb}, {j, 1, OptionValue["Nimp"]}];
			(* Ust n_a,up n_b,dw *)
			If[Ust != 0.0 && Norb > 1,
				Hsector += SparseArray @ DiagonalMatrix[
					Ust * Sum[
						If[orbA != orbB,
							num[[j, 1, orbA]] * num[[j, 2, orbB]],
						(* else *)
							0
						]
					, {orbA, Norb}, {orbB, Norb}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* Usec * (n_a,up n_b,up + n_a,dw n_b,dw) *)
			If[Usec != 0.0 && Norb > 1,
				Hsector += Usec * SparseArray @ DiagonalMatrix[ 
					Sum[
						num[[j, \[Sigma], orb]] * num[[j, \[Sigma], orb+1]]
					, {\[Sigma], f}, {orb, Norb-1}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* - \[Mu]eff * N *)
			If[\[Mu]eff != 0.0,
				Hsector += - \[Mu]eff * SparseArray @ DiagonalMatrix[
					Sum[
						num[[j, \[Sigma], orb]]
					, {\[Sigma], f}, {orb, Norb}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* \[Delta]_a N_a *)
			Do[
				If[\[Delta][[orb]] == 0.0, Continue[]; ];
				Hsector += \[Delta][[orb]] * SparseArray @ DiagonalMatrix[
					Sum[
						num[[j, \[Sigma], orb]]
					, {j, 1, OptionValue["Nimp"]}, {\[Sigma], f}]
				];
			, {orb, Norb}];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"]
];
Options[HLocalNormal] = {"Nimp" -> 1};
(* *)
HLocalSuperc[L_, f_, Norb_, Sectors_, EdMode_, InteractionParameters_, OptionsPattern[]] := Module[{
	U = InteractionParameters[[1]],
	Ust = InteractionParameters[[2]],
	Usec = InteractionParameters[[3]],
	Jph = InteractionParameters[[4]],
	\[CapitalDelta]ext = InteractionParameters[[5]],
	\[Mu]eff = InteractionParameters[[6]],
	\[Delta] = InteractionParameters[[7]],
	H,Hsector,Hblock,rules,dispatch,num,norb,dim,\[Psi]1,\[Chi],rows,cols,pos,\[CapitalSigma]
	},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = SparseArray[{}, {dim, dim}];
			(* exception: if all parameters vanish, avoid calculations *)
			If[# == ConstantArray[0, Length[#]],
				Sow[Hsector, "Hsector"];
				Continue[];,
			(* else, compute n_j,\[Sigma],orb on all the states of the sector *)
				num = Table[
					n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
				, {j, 1, OptionValue["Nimp"]}, {\[Sigma], f}, {orb, Norb}];
			] &@ Flatten[InteractionParameters];
			(* orbital-wise Hubbard term *)
			Do[
				If[U[[orb]] == 0, Continue[]; ];
				norb = Sum[num[[j,\[Sigma],orb]], {\[Sigma],f}];
				Hsector += 0.5*U[[orb]] * SparseArray @ DiagonalMatrix[
					norb * (norb - 1)
				];		
			, {orb, Norb}, {j, 1, OptionValue["Nimp"]}];
			(* Ust n_a,up n_b,dw *)
			If[Ust != 0.0 && Norb > 1,
				Hsector += SparseArray @ DiagonalMatrix[
					Ust * Sum[
						If[orbA != orbB,
							num[[j, 1, orbA]] * num[[j, 2, orbB]],
						(* else *)
							0
						]
					, {orbA, Norb}, {orbB, Norb}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* Usec * (n_a,up n_b,up + n_a,dw n_b,dw) *)
			If[Usec != 0.0 && Norb > 1,
				Hsector += Usec * SparseArray @ DiagonalMatrix[ 
					Sum[
						num[[j, \[Sigma], orb]] * num[[j, \[Sigma], orb+1]]
					, {\[Sigma], f}, {orb, Norb-1}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* pair hopping *)
			Hblock = SparseArray[{}, {dim, dim}];
			Do[
				If[orb2 > orb1,
					If[Jph[[orb1, orb2]] == 0.0, Continue[]; ];
					\[Psi]1 = PairHoppingSelect[L, f, 1, orb1, orb2, \[Psi]];
					If[Length[\[Psi]1] == 0, Continue[];];
					\[Chi] = PairHopping[L, f, 1, orb1, orb2, \[Psi]1];
					rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
					\[CapitalSigma] = CCSign[L, f, {j,j,j,j}, {1,2,2,1}, {orb1,orb1,orb2,orb2}, \[Psi]1];
					Hblock += Jph[[orb1, orb2]] * SparseArray[pos -> \[CapitalSigma], {dim, dim}];
				];
			, {orb1, Norb}, {orb2, Norb}, {j, 1, OptionValue["Nimp"]}];
			Hblock = Hblock + Hblock\[ConjugateTranspose];
			Hsector += Hblock;
			(* external pairing field *)
			Hblock = SparseArray[{}, {dim, dim}];
			Do[
				If[\[CapitalDelta]ext[[orb]] != 0.0,
					\[Psi]1 = CreatePairSelect[L, f, j, j, 1, 2, orb, orb, \[Psi]];
					If[Length[\[Psi]1] != 0, 
						\[Chi] = CreatePair[L, f, j, j, 1, 2, orb, orb, \[Psi]1];
						rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						\[CapitalSigma] = CCSign[L, f, {j,j}, {1,2}, {orb,orb}, \[Psi]1];
						Hblock += \[CapitalDelta]ext[[orb]] * SparseArray[pos -> \[CapitalSigma], {dim, dim}];
					];
				];
			, {orb, Norb}, {j, 1, OptionValue["Nimp"]}];
			Hblock = Hblock + Hblock\[ConjugateTranspose];
			Hsector += Hblock;
			(* - \[Mu]eff * N *)
			If[\[Mu]eff != 0.0,
				Hsector += - \[Mu]eff * SparseArray @ DiagonalMatrix[
					Sum[
						num[[j, \[Sigma], orb]]
					, {\[Sigma], f}, {orb, Norb}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* \[Delta]_a N_a *)
			Do[
				If[\[Delta][[orb]] == 0.0, Continue[]; ];
				Hsector += \[Delta][[orb]] * SparseArray @ DiagonalMatrix[
					Sum[
						num[[j, \[Sigma], orb]]
					, {j, 1, OptionValue["Nimp"]}, {\[Sigma], f}]
				];
			, {orb, Norb}];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"]
];
Options[HLocalSuperc] = {"Nimp" -> 1};
(* *)
HLocalRaman[L_, f_, Norb_, Sectors_, EdMode_, InteractionParameters_, OptionsPattern[]] := Module[{
	\[CapitalOmega] = InteractionParameters[[1]],
	U = InteractionParameters[[2]],
	Ust = InteractionParameters[[3]],
	Usec = InteractionParameters[[4]],
	Jse = InteractionParameters[[5]],
	\[Mu]eff = InteractionParameters[[6]],
	\[Delta] = InteractionParameters[[7]],
	shift = If[OptionValue["HFMode"], 0, (* else *)0.5],
	H,Hsector,Hblock,rules,dispatch,num,norb,dim,\[Psi]1,\[Chi],rows,cols,pos,\[CapitalSigma]
	},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = SparseArray[{}, {dim, dim}];
			(* exception: if all parameters vanish, avoid calculations *)
			If[# == ConstantArray[0, Length[#]],
				Sow[Hsector, "Hsector"];
				Continue[];,
			(* else, compute n_j,\[Sigma],orb on all the states of the sector *)
				num = Table[
					n[L, f, Norb, j, \[Sigma], orb, \[Psi]] - shift
				, {j, 1, OptionValue["Nimp"]}, {\[Sigma], f}, {orb, Norb}];
			] &@ Flatten[InteractionParameters];
			(* Raman term *)
			Do[
				(*If[\[CapitalOmega][[orb, \[Sigma], \[Rho]]] == 0.0, Continue[]; ];*)
				If[
					\[Rho] == \[Sigma], (* diagonal elements of Raman matrix represent a magnetic field *)
					Hsector += \[CapitalOmega][[orb,\[Sigma],\[Sigma]]] * SparseArray @ DiagonalMatrix[num[[j,\[Sigma],orb]] + shift];,	
				(* else, if \[Rho] > \[Sigma], off diagonal elements of Raman matrix *)
					\[Psi]1 = HopSelect[L, f, 1, 1, \[Rho], \[Sigma], orb, orb, \[Psi]];
					If[Length[\[Psi]1] != 0, 
						\[Chi] = Hop[L, f, 1, 1, \[Rho], \[Sigma], orb, orb, \[Psi]1];
						rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						\[CapitalSigma] = CCSign[L, f, {1,1}, {\[Rho],\[Sigma]}, {orb,orb}, \[Psi]1];
						If[Index[L, f, Norb, 1, \[Rho], orb] > Index[L, f, Norb, 1, \[Sigma], orb], 
							\[CapitalSigma] = -\[CapitalSigma]
						];(* if \[Rho] > \[Sigma], we are applying e.g. cdg_2 c_1. When moving cdg_2 before position 2, it jumps over c_1 and the sign changes! *)
						Hblock = \[CapitalOmega][[orb,\[Sigma],\[Rho]]] * SparseArray[pos -> \[CapitalSigma], {dim, dim}];
						Hsector += Hblock + Hblock\[ConjugateTranspose];			
					];
				];
			, {orb, Norb}, {j, 1, OptionValue["Nimp"]}, {\[Rho], f}, {\[Sigma], \[Rho], f}];
			(* orbital-wise Hubbard term *)
			Do[
				If[U[[orb]] == 0, Continue[]; ];
				norb = Sum[num[[j,\[Sigma],orb]], {\[Sigma],f}];
				Hsector += 0.5*U[[orb]] * SparseArray @ DiagonalMatrix[
					If[OptionValue["HFMode"], norb*(norb-1), (* else *)norb^2]
				];		
			, {orb, Norb}, {j, 1, OptionValue["Nimp"]}];
			(* Ust n_a,up n_b,dw *)
			If[Ust != 0.0 && Norb > 1,
				Hsector += SparseArray @ DiagonalMatrix[
					Ust * Sum[
						If[orbA != orbB,
							num[[j, 1, orbA]] * num[[j, 2, orbB]],
						(* else *)
							0
						]
					, {orbA, Norb}, {orbB, Norb}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* Usec * (n_a,up n_b,up + n_a,dw n_b,dw) *)
			If[Usec != 0.0 && Norb > 1,
				Hsector += Usec * SparseArray @ DiagonalMatrix[ 
					Sum[
						num[[j, \[Sigma], orb]] * num[[j, \[Sigma], orb+1]]
					, {\[Sigma], f}, {orb, Norb-1}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* Jse  *)
			If[Jse != 0.0,
				Hblock = SparseArray[{},{dim,dim}];
				Do[
					\[Psi]1 = SpinExchangeSelect[L, f, 1, orb1, orb2, \[Psi]];
					If[Length[\[Psi]1] == 0, Continue[];];
					\[Chi] = SpinExchange[L, f, 1, orb1, orb2, \[Psi]1];
					rows = \[Chi]/.dispatch;(* *)cols = \[Psi]1/.dispatch;(* *)pos = {rows,cols}\[Transpose];
					\[CapitalSigma] = -CCSign[L, f, {j,j,j,j}, {1,2,2,1}, {orb1,orb1,orb2,orb2}, \[Psi]1];
					Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
				, {orb1, Norb}, {orb2, orb1+1, Norb}, {j, 1, OptionValue["Nimp"]}];
				Hsector += Jse * (Hblock + Hblock\[ConjugateTranspose]);
			];
			(* - \[Mu]eff * N *)
			If[\[Mu]eff != 0.0,
				Hsector += - \[Mu]eff * SparseArray @ DiagonalMatrix[
					Sum[
						num[[j, \[Sigma], orb]] + shift
					, {\[Sigma], f}, {orb, Norb}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* \[Delta]_a N_a *)
			Do[
				If[\[Delta][[orb]] == 0.0, Continue[]; ];
				Hsector += \[Delta][[orb]] * SparseArray @ DiagonalMatrix[
					Sum[
						num[[j, \[Sigma], orb]] + shift
					, {j, 1, OptionValue["Nimp"]}, {\[Sigma], f}]
				];
			, {orb, Norb}];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"]
];
Options[HLocalRaman] = {"Nimp" -> 1, "HFMode" -> True};
(* *)
HLocalMagnetic[L_, f_, Norb_, Sectors_, EdMode_, InteractionParameters_, OptionsPattern[]] := Module[{
	\[CapitalOmega] = InteractionParameters[[1]],
	U = InteractionParameters[[2]],
	Ust = InteractionParameters[[3]],
	Usec = InteractionParameters[[4]],
	\[Mu]eff = InteractionParameters[[5]],
	\[Delta] = InteractionParameters[[6]],
	shift = If[OptionValue["HFMode"], 0, (* else *)0.5],
	H,Hsector,Hblock,rules,dispatch,num,norb,dim,\[Psi]1,\[Chi],rows,cols,pos,\[CapitalSigma]
	},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = SparseArray[{}, {dim, dim}];
			(* exception: if all parameters vanish, avoid calculations *)
			If[# == ConstantArray[0, Length[#]],
				Sow[Hsector, "Hsector"];
				Continue[];,
			(* else, compute n_j,\[Sigma],orb on all the states of the sector *)
				num = Table[
					n[L, f, Norb, j, \[Sigma], orb, \[Psi]] - shift
				, {j, 1, OptionValue["Nimp"]}, {\[Sigma], f}, {orb, Norb}];
			] &@ Flatten[InteractionParameters];
			(* Magnetic term *)
			Do[	
				Hsector += \[CapitalOmega][[orb,\[Sigma]]] * SparseArray[ Band[{1,1}] -> num[[j,\[Sigma],orb]] + shift ];	
			, {orb, Norb}, {j, 1, OptionValue["Nimp"]}, {\[Sigma], f}];
			(* orbital-wise Hubbard term *)
			Do[
				If[U[[orb]] == 0, Continue[]; ];
				norb = Sum[num[[j,\[Sigma],orb]], {\[Sigma],f}];
				Hsector += 0.5*U[[orb]] * SparseArray[ 
					Band[{1,1}] -> 
					If[OptionValue["HFMode"], norb*(norb-1), (* else *)norb^2]
				];		
			, {orb, Norb}, {j, 1, OptionValue["Nimp"]}];
			(* Ust n_a,up n_b,dw *)
			If[Ust != 0.0 && Norb > 1,
				Hsector += SparseArray[
					Band[{1,1}] ->
					Ust * Sum[
						If[orbA != orbB,
							num[[j, 1, orbA]] * num[[j, 2, orbB]],
						(* else *)
							0
						]
					, {orbA, Norb}, {orbB, Norb}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* Usec * (n_a,up n_b,up + n_a,dw n_b,dw) *)
			If[Usec != 0.0 && Norb > 1,
				Hsector += Usec * SparseArray[ 
					Band[{1,1}] ->
					Sum[
						num[[j, \[Sigma], orb]] * num[[j, \[Sigma], orb+1]]
					, {\[Sigma], f}, {orb, Norb-1}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* - \[Mu]eff * N *)
			If[\[Mu]eff != 0.0,
				Hsector += - \[Mu]eff * SparseArray[
					Band[{1,1}] ->
					Sum[
						num[[j, \[Sigma], orb]] + shift
					, {\[Sigma], f}, {orb, Norb}, {j, 1, OptionValue["Nimp"]}]
				];
			];
			(* \[Delta]_a N_a *)
			Do[
				If[\[Delta][[orb]] == 0.0, Continue[]; ];
				Hsector += \[Delta][[orb]] * SparseArray[
					Band[{1,1}] ->
					Sum[
						num[[j, \[Sigma], orb]] + shift
					, {j, 1, OptionValue["Nimp"]}, {\[Sigma], f}]
				];
			, {orb, Norb}];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"]
];
Options[HLocalMagnetic] = {"Nimp" -> 1, "HFMode" -> True};

(*  *)
HNonlocal[L_, f_, Norb_, Sectors_, EdMode_] := Which[
	EdMode == "Normal",
	HNonlocalNormal[L, f, Norb, Sectors],
	EdMode == "Superc",
	HNonlocalSuperc[L, f, Norb, Sectors],
	EdMode == "Raman",
	HNonlocalRaman[L, f, Norb, Sectors],
	EdMode == "Magnetic",
	HNonlocalMagnetic[L, f, Norb, Sectors]
];
(*  *)
HLocal[L_, f_, Norb_, Sectors_, EdMode_, InteractionParameters_, OptionsPattern[]] := Module[{
	HFMode = OptionValue["HFMode"],
	Nimp = OptionValue["Nimp"]
	},
	Which[
		EdMode == "Normal",
		HLocalNormal[L, f, Norb, Sectors, EdMode, InteractionParameters],
		EdMode == "Superc",
		HLocalSuperc[L, f, Norb, Sectors, EdMode, InteractionParameters],
		EdMode == "Raman",
		HLocalRaman[L, f, Norb, Sectors, EdMode, InteractionParameters, "Nimp"->Nimp, "HFMode"->HFMode],
		EdMode == "Magnetic",
		HLocalMagnetic[L, f, Norb, Sectors, EdMode, InteractionParameters, "Nimp"->Nimp, "HFMode"->HFMode]
	]
];
Options[HLocal] = {"Nimp" -> 1, "HFMode" -> True}

(* Get the Hamiltonian structure once for all *)
GetHamiltonian[L_, f_, Norb_, Nimp_, HFMode_, Sectors_, LoadHamiltonianQ_, HnonlocFile_, HlocFile_, EdMode_, InteractionParameters_] := Module[
	{HnonlocBlocks, HlocBlocks},
	If[LoadHamiltonianQ,
		Print["Getting Hamiltonians from file"];
		HnonlocBlocks = Import[HnonlocFile];
		HlocBlocks = Import[HlocFile];
		Print["Done! Let's get started! \n"],
	(* else *)
		Print["Computing Hamiltonians..."];
		Print["Time: ", First @ AbsoluteTiming[
			HnonlocBlocks = HNonlocal[L, f, Norb, Sectors, EdMode];
			Export[HnonlocFile, HnonlocBlocks];
			HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode, InteractionParameters, "Nimp"->Nimp, "HFMode"->HFMode];
			Export[HlocFile, HlocBlocks];
		], " sec."];
		Print["Done! Let's get started! \n"];
	];
	{HnonlocBlocks, HlocBlocks}
];

(* build the impurity Hamiltonian from the local and nonlocal blocks and the respective parameters *)
HImp[HnonlocBlocks_, Hloc_, BathParameters_] := Module[
	{FlatBathParameters = Flatten[BathParameters], Hnonloc},	
	(* build non local Hamiltonian -> avoid Dot[] as it returns dense array *)
	Hnonloc = SparseArray[#]&/@(
		Sum[
			FlatBathParameters[[i]]*#[[i]],
		{i, Length@FlatBathParameters}] &/@ HnonlocBlocks
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

(* this is the ladder operator S_\[Sigma]\[Rho] expressed in the basis generated by BuildSector *)
LadderOperator[L_, f_, Norb_, Sectors_, i_, orb_, \[Sigma]_, \[Rho]_] := Module[
	{dim, rules, dispatch, \[Psi]1, \[Psi]2,\[Chi], rows,cols,pos, block},
	Do[
		dim = Length[\[Psi]];
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		block = SparseArray[{},{dim,dim}];
		\[Psi]1 = ChangeSpinSelect[L,f,i,\[Rho],orb,\[Psi]];
		If[
			Length[\[Psi]1] != 0,
			\[Chi]=s[L,f,i,\[Sigma],\[Rho],orb,\[Psi]1];
			rows=\[Chi]/.dispatch;(**)cols=\[Psi]1/.dispatch;(**)pos={rows,cols}\[Transpose];
			block += SparseArray[
				pos -> ConstantArray[1.0, Length[\[Psi]1]]
			,{dim,dim}];
		];
	, {\[Psi],Sectors}];
	block
];


End[]

EndPackage[]
