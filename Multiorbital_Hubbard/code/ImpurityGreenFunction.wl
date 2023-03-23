(* ::Package:: *)

BeginPackage["ImpurityGreenFunction`", {"MyLinearAlgebra`", "Sectors`"}]


GreenFunctionED::usage = "GreenFunctionED[L, f, Norb, {i,j}, \[Sigma], orb, Sectors, QnsSectorList, eigs, T, zlist, EdMode]"

GreenFunctionImpurityNormal::usage = "GreenFunctionImpurity[L_, f_, Norb_, \[Sigma]_, orb_, Egs_, gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_]"

GreenFunctionImpurityNambu::usage = "GreenFunctionImpurityNambu[L_, f_, Norb_, orb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_]"

GreenFunctionImpurityRaman::usage = "GreenFunctionImpurityRaman[L_, f_, Norb_, orb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_]"

GreenFunctionImpurity::usage = "."

InverseGreenFunction::usage = "InverseGreenFunction[L, f, Norb, \[Sigma], orb, Egs, gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist] evaluates numerically the inverse Green function 
of the fully interacting impurity problem. This function returns a list of either numbers when EdMode = Normal, or matrices in the suitable Nambu basis depending on EdMode. The input 
spin and orbital indexes are completely ignored when the Nambu basis mixes up spin / orbital degrees of freedom. "


Begin["Private`"]

Print["Package ImpurityGreenFunction` loaded successfully."];

(*                APPLY CDG / C TO STATES             *)
(* find the quantum number of the final state after application of operator_{j,\[Sigma],orb} to a state living in the sector qns *)
FinalSector[L_, f_, Norb_, j_, \[Sigma]_, orb_, qns_, operator_, EdMode_] := Module[
	{newqns = qns},
	Which[
		operator == "Creation",
		Which[
			EdMode == "Normal",
			If[qns[[f*(orb-1)+\[Sigma]]] == L, Return[]];  (* trivial case: return Null if final sector does not exist *)
			newqns[[f*(orb-1)+\[Sigma]]] += 1;,
		(* ---------------------------- *)
			EdMode == "InterorbNormal",
			If[qns[[\[Sigma]]] == Norb*L, Return[]]; (* trivial case: return Null if sector does not exist *)
			newqns[[\[Sigma]]] += 1;,
		(* ---------------------------- *)
			EdMode == "Raman",
			If[qns[[orb]] == f*L, Return[]]; (* trivial case: return Null if sector does not exist *)
			newqns[[orb]] += 1;,
		(* ---------------------------- *)
			EdMode == "Superc",
			If[f>2, Return["error. f>2 not supported with EdMode = ''Superc''"];];
			If[(qns[[orb]] == -L && \[Sigma] == 2) || (qns[[orb]] == L && \[Sigma] == 1), Return[]];  (* trivial case: return Null if sector does not exist *)
			Which[
				\[Sigma]==1, newqns[[orb]] += 1,
				\[Sigma]==2, newqns[[orb]] -= 1
			];,
		(* ---------------------------- *)
			EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
			If[f>2, Return["error. f>2 not supported with EdMode = ''FullSuperc'' or ''InterorbSuperc''"];];
			If[(qns == -Norb*L && \[Sigma]==2) || (qns == Norb*L && \[Sigma]==1), Return[]];  (* trivial case: return Null if final sector does not exist. *)
			Which[
				\[Sigma]==1, newqns+=1,
				\[Sigma]==2, newqns-=1
			];
		],
	(* ------------------------------------------ *)	
		operator == "Annihilation",
		Which[
			EdMode == "Normal",
			If[qns[[f*(orb-1)+\[Sigma]]] == 0, Return[]];  (* trivial case *)
			newqns[[f*(orb-1)+\[Sigma]]] -= 1;,
		(* ---------------------------- *)
			EdMode == "InterorbNormal",
			If[qns[[\[Sigma]]] == 0, Return[]];  (* trivial case *)
			newqns[[\[Sigma]]] -= 1;,
		(* ---------------------------- *)
			EdMode == "Raman",
			If[qns[[orb]] == 0, Return[]]; (* trivial case: return Null if sector does not exist *)
			newqns[[orb]] -= 1;,
		(* ---------------------------- *)
			EdMode == "Superc",
			If[f>2, Return["error. f>2 not supported with EdMode = ''Superc''"];];
			If[(qns[[orb]] == -L && \[Sigma] == 1) || (qns[[orb]] == L && \[Sigma] == 2), Return[]];  (* trivial case *)
			Which[
				\[Sigma]==1, newqns[[orb]] -= 1,
				\[Sigma]==2, newqns[[orb]] += 1
			];,
		(* ---------------------------- *)
			EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
			If[f>2, Return["error. f>2 not supported with EdMode = ''FullSuperc'' or ''InterorbSuperc''"];];
			If[(qns == -Norb*L && \[Sigma]==1) || (qns == Norb*L && \[Sigma]==2), Return[]];  (* trivial case *)
			Which[
				\[Sigma]==1, newqns -= 1,
				\[Sigma]==2, newqns += 1
			];
		],
	(* ------------------------------------------ *)	
		operator == "Density",
		newqns = qns;	
	];
	newqns
];

(* apply cdg_{j,\[Sigma],orb} |gs>, where |gs> belongs to the sector with quantum numbers qns and give the resulting vector resized to fit the dimension of the sector obtained adding a particle with state label (j,\[Sigma],orb)*)
ApplyCdg[L_, f_, Norb_, j_, \[Sigma]_, orb_, gs_, qns_, Sectors_, SectorsDispatch_, EdMode_] := Module[
	{newqns = qns, startingsector = Sectors[[qns/.SectorsDispatch]],finalsector, newdim,sign,dispatch,pos,newpos,coeff,\[Psi]1,\[Chi]},
	(* build the final sector quantum numbers *)
	Which[
		EdMode == "Normal",
		If[qns[[f*(orb-1)+\[Sigma]]] == L, Return[0]];  (* trivial case *)
		newqns[[f*(orb-1)+\[Sigma]]] += 1;,
	(* ---------------------------- *)
		EdMode == "InterorbNormal",
		If[qns[[\[Sigma]]] == Norb*L, Return[0]];  (* trivial case *)
		newqns[[\[Sigma]]]+=1;,
	(* ---------------------------- *)
		EdMode == "Raman",
		If[qns[[orb]] == f*L, Return[0]]; (* trivial case *)
		newqns[[orb]] += 1;,	
	(* ---------------------------- *)
		EdMode == "Superc",
		If[f>2, Return["error. f>2 not supported with EdMode = ''Superc''"];];
		If[(qns[[orb]] == -L && \[Sigma] == 2) || (qns[[orb]] == L && \[Sigma] == 1),Return[0]];  (* trivial case *)
		Which[
			\[Sigma]==1, newqns[[orb]]+=1,
			\[Sigma]==2, newqns[[orb]]-=1
		],
	(* ---------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		If[f>2, Return["error. f>2 not supported with EdMode = ''FullSuperc'' or ''InterorbSuperc''"];];
		If[(qns==-Norb*L&&\[Sigma]==2)||(qns==Norb*L&&\[Sigma]==1),Return[0]];  (* trivial case *)
		Which[
			\[Sigma]==1, newqns+=1,
			\[Sigma]==2, newqns-=1
		]
	];
	(* check which states of the starting sector can host the extra particle *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,startingsector],1]];
	\[Psi]1 = CreateParticleSelect[L, f, j, \[Sigma], orb, startingsector];
	pos = \[Psi]1/.dispatch;
	(* compute the correct signs obtained moving cdg to the correct position  *)
	sign = CSign[L, f, j, \[Sigma], orb, \[Psi]1];
	(* list of coefficients that remain non vanishing *)
	coeff = gs[[pos]]*sign;
	(* get the final sector *)
	finalsector = Sectors[[newqns/.SectorsDispatch]];
	newdim = Length[finalsector];
	(* create a dispatch that labels all these states *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1]];(* find the new positions where the entries of gs should go after applying cdg_spin *)
	\[Chi] = cdg[L,f,j,\[Sigma],orb,\[Psi]1];
	newpos = \[Chi]/.dispatch;
	(* create a list of rules and define the resulting array *)
	SparseArray[Thread[newpos->coeff],newdim]
];

(* apply c_{j,\[Sigma],orb} |gs>, where |gs> belongs to the sector with quantum numbers qns and give the resulting vector resized to fit the dimension of the sector obtained removing a particle with state label (j,\[Sigma],orb)*)
ApplyC[L_, f_, Norb_, j_, \[Sigma]_, orb_, gs_, qns_, Sectors_, SectorsDispatch_, EdMode_] := Module[
	{newqns = qns,startingsector = Sectors[[qns/.SectorsDispatch]],finalsector, newdim,sign,dispatch,pos,newpos,coeff,\[Psi]1,\[Chi]},
	(* build the final sector quantum numbers *)
	Which[
		EdMode == "Normal",
		If[qns[[f*(orb-1)+\[Sigma]]]==0, Return[0]];  (* trivial case *)
		newqns[[f*(orb-1)+\[Sigma]]]-=1;,
	(* ---------------------------- *)
		EdMode == "InterorbNormal",
		If[qns[[\[Sigma]]]==0, Return[0]];  (* trivial case *)
		newqns[[\[Sigma]]]-=1;,
	(* ---------------------------- *)
		EdMode == "Raman",
		If[qns[[orb]] == 0, Return[0]]; (* trivial case *)
		newqns[[orb]] -= 1;,
	(* ---------------------------- *)
		EdMode == "Superc",
		If[f>2, Return["error. f>2 not supported with EdMode = ''Superc''"];];
		If[(qns[[orb]]==-L&&\[Sigma]==1)||(qns[[orb]]==L&&\[Sigma]==2), Return[0]];  (* trivial case *)
		Which[
			\[Sigma]==1, newqns[[orb]]-=1,
			\[Sigma]==2, newqns[[orb]]+=1
		],
	(* ---------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		If[f>2, Return["error. f>2 not supported with EdMode = ''FullSuperc'' or ''InterorbSuperc''"];];
		If[(qns==-Norb*L&&\[Sigma]==1)||(qns==Norb*L&&\[Sigma]==2),Return[0]];  (* trivial case *)
		Which[
			\[Sigma]==1, newqns-=1,
			\[Sigma]==2, newqns+=1
		]
	];
	(* check which states of the starting sector can host the extra particle *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,startingsector],1]];
	\[Psi]1 = DestroyParticleSelect[L, f, j, \[Sigma], orb, startingsector];
	pos = \[Psi]1/.dispatch;
	(* compute the correct signs obtained moving cdg to the correct position  *)
	sign = CSign[L,f,j,\[Sigma],orb,\[Psi]1];
	(* list of coefficients that remain non vanishing *)
	coeff = gs[[pos]]*sign;
	(* find final sector *)
	finalsector = Sectors[[newqns/.SectorsDispatch]];
	newdim = Length[finalsector];
	(* create a dispatch that labels all these states *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1]];(* find the new positions where the entries of gs should go after applying cdg_spin *)
	\[Chi] = c[L,f,j,\[Sigma],orb,\[Psi]1];
	newpos = \[Chi]/.dispatch;
	(* create a list of rules and define the resulting array *)
	SparseArray[Thread[newpos->coeff], newdim]
];



(*      EXACT CALCULATION OF GREEN FUNCTION       *)
(* compute the Green function using the full spectrum of the given Hamiltonian *)
GreenFunctionED[L_, f_, Norb_, {i_,j_}, {\[Sigma]1_, \[Sigma]2_}, {orb1_, orb2_}, Sectors_, QnsSectorList_, eigs_, T_, zlist_, EdMode_, OptionsPattern[]] := Module[
	{newqns, SectorsDispatch, startingsector, finalsector, dim, newdim, dispatch, sign, rows, cols, pos, \[Psi]1, \[Psi]2, \[Chi]1, \[Chi]2, A, B, P, Q, S, Z=1., e0, energies, eigenstates, G},
	(* get energies and eigenstates *)
	energies = eigs[[All, 1]];
	eigenstates = eigs[[All, 2]];
	e0 = Min[Flatten[energies]];
	(* get normalization factor if required *)
	If[T != 0, 
		Z = Total[ Exp[-(Flatten[energies]-e0)/T] ]
	]; 
	(* get dispatch of the quantum numbers and initialize G.F. *)
	SectorsDispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,QnsSectorList],1]];
	G = ConstantArray[0, Length[zlist]];
	(* loop over all the sectors *)
	Do[
		(* get starting sector *)
		startingsector = Sectors[[qns/.SectorsDispatch]];
		dim = Length[startingsector];
		dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&, startingsector],1]];
		(* select those states for which you can create a particle with index j,\[Sigma]2,orb2 *)
		\[Psi]2 = CreateParticleSelect[L, f, j, \[Sigma]2, orb2, startingsector];
		If[Length[\[Psi]2] == 0, Continue[]; ]; (* if there are some final states *)
			cols = \[Psi]2/.dispatch;
			sign = CSign[L, f, j, \[Sigma]2, orb2, \[Psi]2];
			(* get final sector *)
			newqns = FinalSector[L, f, Norb, 1, \[Sigma]2, orb2, qns, "Creation", EdMode];
			If[newqns == Null, Continue[];];
			finalsector = Sectors[[newqns/.SectorsDispatch]];
			newdim = Length[finalsector];
			dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1]];
			\[Chi]2 = cdg[L,f,j,\[Sigma]2,orb2,\[Psi]2];
			rows = \[Chi]2/.dispatch;
			(* build the matrix Smn = <n|c_i,\[Sigma]1,orb1|m><m|cdg_j,\[Sigma]2,orb2|n> * (e^-\[Beta]Em + e^-\[Beta]En)/(z + En - Em) *)
			(* where |n> is the n-th eigenstate that belongs to the sector with quantum number qns, and |m> is in a sector with one more particle. *)
			(* let Amn = <m|cdg_j,\[Sigma]2,orb2|n> and Bnm = <n|c_i,\[Sigma]1,orb1|m> *)
			pos = {rows, cols}\[Transpose];
			A = SparseArray[Thread[pos->sign], {newdim, dim}];(* A in the canonical basis *)
			P = Transpose[eigenstates[[qns/.SectorsDispatch]]];(* change of basis in starting sector *)
			Q = Transpose[eigenstates[[newqns/.SectorsDispatch]]];(* change of basis in final sector *)
			A = Q\[HermitianConjugate] . A . P;(* A in the eigenstates basis *)
	(* ------------------------------------- *)
		\[Psi]1 = DestroyParticleSelect[L, f, i, \[Sigma]1, orb1, finalsector];
		If[\[Psi]1 == 0, Continue[]; ];
			rows = \[Psi]1/.dispatch;
			sign = CSign[L,f,i,\[Sigma]1,orb1,\[Psi]1];
			dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,startingsector],1]];(* find the new positions where the entries of gs should go after applying cdg_spin *)
			\[Chi]1 = c[L,f,i,\[Sigma]1,orb1,\[Psi]1];
			cols = \[Chi]1/.dispatch;
			pos = {rows, cols}\[Transpose];
			B = SparseArray[Thread[pos->sign], {newdim, dim}];(* B in the canonical basis *)
			B = Q\[HermitianConjugate] . B . P;
	(* ------------------------------------- *)
		S = A * B *
			Table[
				If[T == 0, 
					KroneckerDelta[energies[[qns/.SectorsDispatch]][[n]], e0] + KroneckerDelta[energies[[newqns/.SectorsDispatch]][[m]], e0],
				(* else, if T != 0 *)
					Exp[-(energies[[qns/.SectorsDispatch]][[n]]-e0)/T] + Exp[-(energies[[newqns/.SectorsDispatch]][[m]]-e0)/T]
				]
			, {m, newdim}, {n, dim}];
		G = G + (Total[#,2]&/@(
			S * Table[
				1./(# + energies[[qns/.SectorsDispatch]][[n]] - energies[[newqns/.SectorsDispatch]][[m]])
			, {m, newdim}, {n, dim}] &/@ zlist));
	, {qns, QnsSectorList}];
	G/Z
];
Options[GreenFunctionED] = {NormalizedFunction -> False};



(*
GreenFunctionImpurityNambuOld[L_, f_, Norb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_, OptionsPattern[] ] := Module[{
    NMatsubara = Length[zlist], zlistextended = Join[zlist, -zlist], orb = OptionValue[Orb], OrbitalSymmetry = OptionValue[OrbitalSymmetry],
    GF, GFO, GFP
    },
    Which[ 
        EdMode == "Superc",
        (* initialize Green function as a NMatsubara x 2 x 2 tensor *)
        GF = ConstantArray[0, {NMatsubara, 2, 2}];
        (* if there is spin symmetry, then perform Lanczos only once, apply to +z and -z and then split the result *)
        {GF[[All,1,1]], GF[[All,2,2]]} = Partition[ 
            GreenFunctionImpurity[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
        NMatsubara]; 
        GF[[All,2,2]] = -GF[[All,2,2]]; (* GF_11 -> G_{up,orb1; up,orb1}(z) ;  GF_22 -> - G_{dw,orb1; dw,orb1}(-z) *) 
        (* compute GFO, where Odg = adg_up + a_dw *)
        GFO = GreenFunctionImpurity[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
        (* compute GFP, where Pdg = adg_up + I a_dw *)
        GFP = GreenFunctionImpurity[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
        (* compute the off-diagonal part using the diagonal part and GFO, GFP *)
	    GF[[All,1,2]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All,1,1]] + GF[[All,2,2]]));
	    GF[[All,2,1]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All,1,1]] + GF[[All,2,2]]));,
    (* ------------------------------------------------------------- *)
        EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
        (* initialize Green function as a NMatsubara x 2 x 2 tensor *)
        GF = ConstantArray[0, {NMatsubara, 2*Norb, 2*Norb}];
        If[
        (* general case: NO orbital symmetry *)
            OrbitalSymmetry,
            (* compute the diagonal part of all the Nambu blocks *)
            Do[{
                GF[[All, 2(orb1-1)+1, 2(orb2-1)+1]], (* element 11, 13, 15, 33, 35, 55, ... *)
                GF[[All, 2(orb1-1)+2, 2(orb2-1)+2]] (* element 22, 24, 26, 44, 46, 66, ... *)
                } = Partition[ 
                    GreenFunctionImpurity[L, f, Norb, {orb1,orb2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
                NMatsubara]; 
                GF[[All, 2(orb1-1)+2, 2(orb2-1)+2]] = -GF[[All, 2(orb1-1)+2, 2(orb2-1)+2]]; (* put correct sign on elements 22, 24, 44, ...*)
            , {orb1, Norb}, {orb2, orb1, Norb}];
            Do[
				(* compute GFO, where Odg = cdg_{up,orb1} + c_{dw,orb2} *)
				GFO = GreenFunctionImpurity[L, f, Norb, {orb1,orb2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
				(* compute GFP, where Pdg = cdg_{up,orb1} + I c_{dw,orb2} *)
				GFP = GreenFunctionImpurity[L, f, Norb, {orb1,orb2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
				(* compute the top-right element of each block using the diagonal part and GFO, GFP *)
				GF[[All, 2(orb1-1)+1, 2(orb2-1)+2]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 2(orb1-1)+1, 2(orb1-1)+1]] + GF[[All, 2(orb2-1)+2, 2(orb2-1)+2]]));
				(* compute the bottom-left element of each non-diagonal block *)
				If[orb2 > orb1, 
					GFO = GreenFunctionImpurity[L, f, Norb, {orb2,orb1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
					GFP = GreenFunctionImpurity[L, f, Norb, {orb2,orb1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
					GF[[All, 2(orb1-1)+2, 2(orb2-1)+1]] = Conjugate[
                        (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 2(orb2-1)+1, 2(orb2-1)+1]] + GF[[All, 2(orb1-1)+2, 2(orb1-1)+2]]))
                    ]; (* you have to take the conjugate! For example GF23 is not Fba, but Fba* !  *)
				];
            , {orb1, Norb}, {orb2, orb1, Norb}];
            (* fill up the lower triangular part by conjugating the upper triangular part *)
            Do[
                GF[[All, i, j]] = Conjugate[GF[[All, j, i]]];
            , {i, 2Norb}, {j, i-1}];,
        (* ------------------------------------- *)
        (* else *)
            !OrbitalSymmetry,
            (* compute the diagonal part of (a representative of) the diagonal Nambu block *)
            {GF[[All, 1, 1]], GF[[All, 2, 2]] } = 
            Partition[ 
                GreenFunctionImpurity[L, f, Norb, {1,1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
            NMatsubara]; 
            (* possibly GF13 = GF24 = 0 ? For now let's be safe ... *)
            (* compute the diagonal part of (a representative of) the diagonal Nambu block *)
            {GF[[All, 1, 3]], GF[[All, 2, 4]] } = 
            0*Partition[ 
                GreenFunctionImpurity[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
            NMatsubara]; 
            GF[[All, 2, 4]] = -GF[[All, 2, 4]];
            (* *)
            (* compute GFO, where Odg = cdg_{up,1} + c_{dw,2} *)
			GFO = GreenFunctionImpurity[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
			(* compute GFP, where Pdg = cdg_{up,1} + I c_{dw,2} *)
			GFP = GreenFunctionImpurity[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
            (* compute the top-right element of each block using the diagonal part and GFO, GFP *)
			GF[[All, 1, 2]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 1, 1]] + GF[[All, 2, 2]]));
            GF[[All, 1, 4]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 1, 1]] + GF[[All, 2, 2]]));
            (* possibly GF23 = - GF14 ? For now let's be safe ... *)
            (* compute GFO, where Odg = cdg_{up,1} + c_{dw,2} *)
			GFO = GreenFunctionImpurity[L, f, Norb, {2,1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
			(* compute GFP, where Pdg = cdg_{up,1} + I c_{dw,2} *)
			GFP = GreenFunctionImpurity[L, f, Norb, {2,1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
            GF[[All, 2, 3]] = Conjugate[
                (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 1,1]] + GF[[All, 2, 2]]))
            ];
            (* fill up all the other entries according to the symmetries *)
            (*
            Do[
                GF[[All, 2(orb-1)+1, 2(orb-1)+1]] = GF[[All, 1, 1]];
                GF[[All, 2(orb-1)+2, 2(orb-1)+2]] = GF[[All, 2, 2]];
            , {orb, 2, Norb}]; (* diagonal elements of diagonal blocks *)
            Do[
                If[orb1 == 1 && orb2 == 2, Continue[]; (* this is done already! *) ];
                GF[[All, 2(orb1-1)+1, 2(orb2-1)+1]] = GF[[All, 1, 3]];
                GF[[All, 2(orb1-1)+2, 2(orb2-1)+2]] = GF[[All, 2, 4]];
                GF[[All, 2(orb1-1)+1, 2(orb2-1)+2]] = GF[[All, 1, 4]];
                GF[[All, 2(orb1-1)+2, 2(orb2-1)+3]] = GF[[All, 2, 3]];
            , {orb1, Norb}, {orb2, orb1+1, Norb}]; (* elements of non-diagonal blocks *)
            Do[
                GF[[All, i, j]] = Conjugate[GF[[All, j, i]]];
            , {i, 2Norb}, {j, i-1}]; (* fill up the lower triangular part by conjugating the upper triangular part *)
            *)
        ]
    ];
    GF
];
Options[GreenFunctionImpurityNambu] = {Orb -> 1, OrbitalSymmetry -> False}; (* if orbitals are not symmetric and EdMode = "Superc", you can specify the orbital index *)
*)


(*     LANCZOS CALCULATION OF THE GREEN FUNCTION         *)
(* compute the Green function < O 1/(z-H) Odg > + < Odg 1/(z+H) O > *)
(* where O = c1 c_{orb1,up} + c2 cdg_{orb1,dw} + c3 c_{orb2,up} + c4 cdg_{orb2,dw} *)
GreenFunctionImpurityNormal[L_, f_, Norb_, {orb1_,orb2_}, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_, OptionsPattern[] ] := Module[{
	(* compute cdg_{s,orb}|gs> and c_{s,orb}|gs> for s = up, dw *)
	adgup, bdw, aup, bdgdw, bdgup, bup, adgdw, adw,
	c1 = OptionValue[c1], c2 = OptionValue[c2], c3 = OptionValue[c3], c4 = OptionValue[c4],
	Odggs, Ogs, newqns, H, E0, a, b, GFOparticle, GFOhole, GFO},	
	If[c1 != 0,
		adgup = ApplyCdg[L, f, Norb, 1, 1, orb1, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
		aup = ApplyC[L, f, Norb, 1, 1, orb1, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
	];
	If[c2 != 0,
		adw = ApplyC[L, f, Norb, 1, 2, orb1, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
		adgdw = ApplyCdg[L, f, Norb, 1, 2, orb1, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
	];
	If[c3 != 0,
		bdgup = ApplyCdg[L, f, Norb, 1, 1, orb2, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
		bup = ApplyC[L, f, Norb, 1, 1, orb2, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
	];
	If[c4 != 0,
		bdw = ApplyC[L, f, Norb, 1, 2, orb2, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
		bdgdw = ApplyCdg[L, f, Norb, 1, 2, orb2, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
	];
(* compute all the main contributions to the Green function *)
(*          G_O(z) "Particle" contribution             *)
	Odggs = c1*adgup + c2*adw + c3*bdgup + c4*bdw; 
	newqns = FinalSector[L, f, Norb, 1, 1, orb1, GsQns, "Creation", EdMode]; (* evaluate the quantum numbers of the final sector *)
	If[newqns === Null, (* if there is no final state, this does not contribute to the GF *)
		GFOparticle = ConstantArray[0.0+0.0*I, Length[zlist]];, 
	(* else *)
		H = Hsectors[[newqns/.SectorsDispatch]]; (* Hamiltonian on that sector *)
		If[Length[H] == 1, (* if the Hamiltonian in the final sector is just a number, avoid Lanczos *)
			GFOparticle = ((Norm[Odggs]^2)/(# - H[[1,1]] + Egs)) &/@ zlist;,
		(* else *)
			{E0,a,b} = Lanczos[H, Normalize[Odggs] ]; (* Apply Lanczos starting from Odg|gs> *)
			H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a] ]; (* Krylov matrix in the final sector *)
			GFOparticle = (Norm[Odggs]^2)*(
				InverseElement[
					SparseArray[(# + Egs) * IdentityMatrix[Length[a] ] - H]
				, {1, 1}] &/@ zlist);
		];
	];
(*           G_O(z) "Hole" contribution               *)
	Ogs = (c1\[Conjugate])*aup + (c2\[Conjugate])*adgdw + (c3\[Conjugate])*bup + (c4\[Conjugate])*bdgdw; 
	newqns = FinalSector[L, f, Norb, 1, 1, orb1, GsQns, "Annihilation", EdMode]; (* evaluate the quantum numbers of the final sector *)
	If[newqns === Null, (* if there is no final state, this does not contribute to the GF *)
		GFOhole = ConstantArray[0.0+0.0*I, Length[zlist]];, 
	(* else *)
		H = Hsectors[[newqns/.SectorsDispatch]]; (* Hamiltonian on that sector *)
		If[Length[H] == 1, (* if the Hamiltonian in the final sector is just a number, avoid Lanczos *)
			GFOhole = ((Norm[Ogs]^2)/(# + H[[1,1]] - Egs)) &/@ zlist;,
		(* else *)
			{E0,a,b} = Lanczos[H, Normalize[Ogs] ]; (* Apply Lanczos starting from O|gs> *)
			H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a] ]; (* Krylov matrix in the final sector *)
			GFOhole = (Norm[Ogs]^2)*(
				InverseElement[
					SparseArray[(# - Egs) * IdentityMatrix[Length[a] ] + H]
				, {1, 1}] &/@ zlist);
		];
	];
	GFOparticle + GFOhole
];
Options[GreenFunctionImpurityNormal] = {c1 -> 1.0, c2 -> 0.0, c3 -> 0.0, c4 -> 0.0}; (* by default use Odg = adgup *)



(* compute the Green function in the Nambu formalism *)
GreenFunctionImpurityNambu[L_, f_, Norb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_, OptionsPattern[] ] := Module[{
    NMatsubara = Length[zlist], zlistextended = Join[zlist, -zlist], orb = OptionValue[Orb], OrbitalSymmetry = OptionValue[OrbitalSymmetry],
    GF, GFO, GFP
    },
    Which[ 
        EdMode == "Superc",
        (* initialize Green function as a NMatsubara x 2 x 2 tensor *)
        GF = ConstantArray[0, {NMatsubara, 2, 2}];
        (* if there is spin symmetry, then perform Lanczos only once, apply to +z and -z and then split the result *)
        {GF[[All,1,1]], GF[[All,2,2]]} = Partition[ 
            GreenFunctionImpurityNormal[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
        NMatsubara]; 
        GF[[All,2,2]] = -GF[[All,2,2]]; (* GF_11 -> G_{up,orb1; up,orb1}(z) ;  GF_22 -> - G_{dw,orb1; dw,orb1}(-z) *) 
        (* compute GFO, where Odg = adg_up + a_dw *)
        GFO = GreenFunctionImpurityNormal[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.0];
        (* compute GFP, where Pdg = adg_up + I a_dw *)
        GFP = GreenFunctionImpurityNormal[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.0*I];
        (* compute the off-diagonal part using the diagonal part and GFO, GFP *)
	    GF[[All,1,2]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All,1,1]] + GF[[All,2,2]]));
	    GF[[All,2,1]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All,1,1]] + GF[[All,2,2]]));,
    (* ------------------------------------------------------------- *)
        EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
        (* initialize Green function as a NMatsubara x 2 x 2 tensor *)
        GF = ConstantArray[0, {NMatsubara, 2*Norb, 2*Norb}];
        (* compute the diagonal part of (a representative of) the diagonal Nambu block *)
        {GF[[All, 1, 1]], GF[[All, 2, 2]] } = 
        Partition[ 
            GreenFunctionImpurityNormal[L, f, Norb, {1,1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
        NMatsubara];    
        GF[[All, 2, 2]] = -GF[[All, 2, 2]];
        If[OrbitalSymmetry,
            GF[[All, 3, 3]] = GF[[All, 1, 1]];
            GF[[All, 4, 4]] = GF[[All, 2, 2]];,
        (* else *)
            {GF[[All, 3, 3]], GF[[All, 4, 4]] } = 
            Partition[ 
                GreenFunctionImpurityNormal[L, f, Norb, {2,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
            NMatsubara];
            GF[[All, 4, 4]] = -GF[[All, 4, 4]];
        ];
        (* 1,2 *)
        GFO = GreenFunctionImpurityNormal[L, f, Norb, {1,1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.0];
        GFP = GreenFunctionImpurityNormal[L, f, Norb, {1,1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.0*I];
	    GF[[All,1,2]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All,1,1]] + GF[[All,2,2]]));
	    GF[[All,2,1]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All,1,1]] + GF[[All,2,2]]));
	    (* 3,4 *)
	    If[OrbitalSymmetry,
	        GF[[All, 3, 4]] = GF[[All, 1, 2]];
	        GF[[All, 4, 3]] = GF[[All, 2, 1]];,
	    (* else *)
            GFO = GreenFunctionImpurityNormal[L, f, Norb, {2,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c1 -> 0., c3 -> 1.0, c4 -> 1.0];
            GFP = GreenFunctionImpurityNormal[L, f, Norb, {2,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c1 -> 0., c3 -> 1.0, c4 -> 1.0*I];
	        GF[[All, 3, 4]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 3, 3]] + GF[[All, 4, 4]]));
	        GF[[All, 4, 3]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All, 3, 3]] + GF[[All, 4, 4]]));
	    ];
	    (* 1,3 *)
        GFO = GreenFunctionImpurityNormal[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c3 -> 1.0];
        GFP = GreenFunctionImpurityNormal[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c3 -> 1.0*I];
	    GF[[All,1,3]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All,1,1]] + GF[[All,3,3]]));
	    GF[[All,3,1]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All,1,1]] + GF[[All,3,3]]));
	    (* 2,4 *)
	    GFO = GreenFunctionImpurityNormal[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c1 -> 0.0, c2 -> 1.0, c4 -> 1.0];
        GFP = GreenFunctionImpurityNormal[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c1 -> 0.0, c2 -> 1.0, c4 -> 1.0*I];
	    GF[[All,2,4]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All,2,2]] + GF[[All,4,4]]));
	    GF[[All,4,2]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All,2,2]] + GF[[All,4,4]]));
        (* 1,4 *)
        GFO = GreenFunctionImpurityNormal[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c1 -> 1.0, c4 -> 1.0];
        GFP = GreenFunctionImpurityNormal[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c1 -> 1.0, c4 -> 1.0*I];
	    GF[[All,1,4]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All,1,1]] + GF[[All,4,4]]));
	    GF[[All,4,1]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All,1,1]] + GF[[All,4,4]]));
	    (* 2,3 *)
        GFO = GreenFunctionImpurityNormal[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c1 -> 0.0, c2 -> 1.0, c3 -> 1.0];
        GFP = GreenFunctionImpurityNormal[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c1 -> 0.0, c2 -> 1.0, c3 -> 1.0*I];
	    GF[[All,2,3]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All,2,2]] + GF[[All,3,3]]));
	    GF[[All,3,2]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All,2,2]] + GF[[All,3,3]]));
    ];
    GF
];
Options[GreenFunctionImpurityNambu] = {Orb -> 1, OrbitalSymmetry -> False};


(* compute the Green function using Raman formalism *)
(* compute the Green function < c_{orb, \[Sigma]} 1/(z-H) cdg_{orb, \[Sigma]} > + < cdg_{orb, \[Sigma]} 1/(z+H) c_{orb, \[Sigma]} > *)
(* where O = c1 c_{orb, \[Sigma]} + c2 c_{orb, \[Rho]} *)
GreenFunctionImpurityNormalRaman[L_, f_, Norb_, orb_, {\[Sigma]_, \[Rho]_}, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_, OptionsPattern[]] := Module[
	{
	(* compute cdg_{orb, \[Sigma]}|gs> and c_{orb, \[Rho]}|gs> *)
	adg\[Sigma], a\[Sigma], adg\[Rho], a\[Rho],
	c1 = OptionValue[c1], c2 = OptionValue[c2],
	Odggs, Ogs, newqns, H, E0, a, b, GFOparticle, GFOhole, GFO
	},	
	If[c1 != 0,
		adg\[Sigma] = ApplyCdg[L, f, Norb, 1, \[Sigma], orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
		a\[Sigma] = ApplyC[L, f, Norb, 1, \[Sigma], orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
	];
	If[c2 != 0,
		adg\[Rho] = ApplyCdg[L, f, Norb, 1, \[Rho], orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
		a\[Rho] = ApplyC[L, f, Norb, 1, \[Rho], orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
	];
(* compute all the main contributions to the Green function *)
(*          G_O(z) "Particle" contribution             *)
	Odggs = c1 * adg\[Sigma] + c2 * adg\[Rho];
	newqns = FinalSector[L, f, Norb, 1, \[Sigma], orb, GsQns, "Creation", EdMode]; (* evaluate the quantum numbers of the final sector *)
	If[newqns === Null, (* if there is no final state, this does not contribute to the GF *)
		GFOparticle = ConstantArray[0.0+0.0*I, Length[zlist]];, 
	(* else *)
		H = Hsectors[[newqns/.SectorsDispatch]]; (* Hamiltonian on that sector *)
		If[Length[H] == 1, (* if the Hamiltonian in the final sector is just a number, avoid Lanczos *)
			GFOparticle = ((Norm[Odggs]^2)/(# - H[[1,1]] + Egs)) &/@ zlist;,
		(* else *)
			{E0,a,b} = Lanczos[H, Normalize[Odggs] ]; (* Apply Lanczos starting from Odg|gs> *)
			H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a] ]; (* Krylov matrix in the final sector *)
			GFOparticle = (Norm[Odggs]^2)*(
				InverseElement[
					SparseArray[(# + Egs) * IdentityMatrix[Length[a] ] - H]
				, {1, 1}] &/@ zlist);
		];
	];
(*           G_O(z) "Hole" contribution               *)
	Ogs = (c1\[Conjugate])*a\[Sigma] + (c2\[Conjugate])*a\[Rho]; 
	newqns = FinalSector[L, f, Norb, 1, \[Sigma], orb, GsQns, "Annihilation", EdMode]; (* evaluate the quantum numbers of the final sector *)
	If[newqns === Null, (* if there is no final state, this does not contribute to the GF *)
		GFOhole = ConstantArray[0.0+0.0*I, Length[zlist]];, 
	(* else *)
		H = Hsectors[[newqns/.SectorsDispatch]]; (* Hamiltonian on that sector *)
		If[Length[H] == 1, (* if the Hamiltonian in the final sector is just a number, avoid Lanczos *)
			GFOhole = ((Norm[Ogs]^2)/(# + H[[1,1]] - Egs)) &/@ zlist;,
		(* else *)
			{E0,a,b} = Lanczos[H, Normalize[Ogs] ]; (* Apply Lanczos starting from O|gs> *)
			H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a] ]; (* Krylov matrix in the final sector *)
			GFOhole = (Norm[Ogs]^2)*(
				InverseElement[
					SparseArray[(# - Egs) * IdentityMatrix[Length[a] ] + H]
				, {1, 1}] &/@ zlist);
		];
	];
	GFOparticle + GFOhole
];
Options[GreenFunctionImpurityNormalRaman] = {c1 -> 1.0, c2 -> 0.0}; (* by default Odg = cdg_{orb,\[Sigma]} *)

(* compute the Green function in the Raman formalism *)
GreenFunctionImpurityRaman[L_, f_, Norb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_, OptionsPattern[] ] := Module[{
	NMatsubara = Length[zlist], orb = OptionValue[Orb], OrbitalSymmetry = OptionValue[OrbitalSymmetry], GF, GFO, GFP
    },
	(* initialize Green function as a NMatsubara x f x f tensor *)
    GF = ConstantArray[0, {NMatsubara, f, f}];
    (* compute diagonal component of the tensor *)
    Table[
         GF[[All, \[Sigma], \[Sigma]]] = GreenFunctionImpurityNormalRaman[L, f, Norb, orb, {\[Sigma], \[Sigma]}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist];
    , {\[Sigma], f}];
    (* compute off-diagonal component of the tensor *)
    Table[
         (* compute GFO, where Odg = adg_\[Sigma] + adg_\[Rho] *)
         GFO = GreenFunctionImpurityNormalRaman[L, f, Norb, orb, {\[Sigma], \[Rho]}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.0];
         (* compute GFP, where Pdg = adg_\[Sigma] + I adg_\[Rho] *)
         GFP = GreenFunctionImpurityNormalRaman[L, f, Norb, orb, {\[Sigma], \[Rho]}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> -1.0];
         (* compute the off-diagonal part using the diagonal part and GFO, GFP *)
	     (*GF[[All, \[Rho], \[Sigma]]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, \[Sigma], \[Sigma]]] + GF[[All, \[Rho], \[Rho]]]));
	     GF[[All, \[Sigma], \[Rho]]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All, \[Sigma], \[Sigma]]] + GF[[All, \[Rho], \[Rho]]]));*)
	     GF[[All, \[Rho], \[Sigma]]] = (1./2.) * (GFO - GF[[All,\[Sigma],\[Sigma]]] - GF[[All,\[Rho],\[Rho]]]);
	     GF[[All, \[Sigma], \[Rho]]] = GF[[All,\[Rho],\[Sigma]]];
    , {\[Sigma], f}, {\[Rho], \[Sigma]+1, f}];
    GF
];
Options[GreenFunctionImpurityRaman] = {Orb -> 1, OrbitalSymmetry -> True};

(*
GreenFunctionImpurityRaman[L_, f_, Norb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_, OptionsPattern[] ] := Module[
	{adggs, ags, newqns, H, E0, a, b, orb = OptionValue[Orb], GFOparticle, GFOhole, GF = ConstantArray[0.+0.*I, {Length[zlist], f, f}]},	
	Do[
		adggs = ApplyCdg[L, f, Norb, 1, \[Rho], orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
		ags = ApplyC[L, f, Norb, 1, \[Sigma], orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode];
		Print[adggs];
		Print[ags];
(*           G_O(z) "Particle" contribution               *)
		newqns = FinalSector[L, f, Norb, 1, \[Rho], orb, GsQns, "Creation", EdMode]; (* evaluate the quantum numbers of the final sector *)
		If[newqns === Null, (* if there is no final state, this does not contribute to the GF *)
			GFOparticle = ConstantArray[0.0+0.0*I, Length[zlist]];, 
		(* else *)
			H = Hsectors[[newqns/.SectorsDispatch]]; (* Hamiltonian on that sector *)
			If[Length[H] == 1, (* if the Hamiltonian in the final sector is just a number, avoid Lanczos *)
				GFOparticle = ((Norm[adggs]^2)/(# - H[[1,1]] + Egs)) &/@ zlist;,
			(* else *)
				{E0,a,b} = Lanczos[H, Normalize[adggs] ]; (* Apply Lanczos starting from Odg|gs> *)
				H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a] ]; (* Krylov matrix in the final sector *)
				GFOparticle = (Norm[adggs]^2)*(
					InverseElement[
						SparseArray[(# + Egs) * IdentityMatrix[Length[a] ] - H]
					, {1, 1}] &/@ zlist);
			];
		];
(*           G_O(z) "Hole" contribution               *) 
		newqns = FinalSector[L, f, Norb, 1, \[Sigma], orb, GsQns, "Annihilation", EdMode]; (* evaluate the quantum numbers of the final sector *)
		If[newqns === Null, (* if there is no final state, this does not contribute to the GF *)
			GFOhole = ConstantArray[0.0+0.0*I, Length[zlist]];, 
		(* else *)
			H = Hsectors[[newqns/.SectorsDispatch]]; (* Hamiltonian on that sector *)
			If[Length[H] == 1, (* if the Hamiltonian in the final sector is just a number, avoid Lanczos *)
				GFOhole = ((Norm[ags]^2)/(# + H[[1,1]] - Egs)) &/@ zlist;,
			(* else *)
				{E0,a,b} = Lanczos[H, Normalize[ags] ]; (* Apply Lanczos starting from O|gs> *)
				H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a] ]; (* Krylov matrix in the final sector *)
				GFOhole = (Norm[ags]^2)*(
					InverseElement[
						SparseArray[(# - Egs) * IdentityMatrix[Length[a] ] + H]
					, {1, 1}] &/@ zlist);
			];
		];
		GF[[All, \[Sigma], \[Rho]]] = GFOparticle + GFOhole;
	, {\[Sigma], f}, {\[Rho], f}];
	GF
];
Options[GreenFunctionImpurityRaman] = {Orb -> 1};
*)


(* Evaluate impurity Green function calling the right function depending on EdMode *)
GreenFunctionImpurity[L_, f_, Norb_, \[Sigma]_, orb_, Egs_, gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_] := 
	Which[
		EdMode == "Normal", 
		GreenFunctionImpurityNormal[L, f, Norb, {orb,orb}, Egs, gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist],
	(* ------------------------------------------- *)
		EdMode == "Raman",
		GreenFunctionImpurityRaman[L, f, Norb, Egs, gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, Orb -> orb],
	(* ------------------------------------------- *)
		EdMode == "Superc" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc", 
		GreenFunctionImpurityNambu[L, f, Norb, Egs, gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, Orb -> orb]
	];

(* invert the Green function depending on EdMode *)
InverseGreenFunction[G_, EdMode_] := 
	Which[
		EdMode == "Normal", 1./G,
		EdMode == "Raman", If[Length[G[[1]]] == 2, TwoByTwoInverse[G], (* else *) Inverse /@ G],
		EdMode == "Superc", TwoByTwoInverse[G],
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc", Inverse /@ G
	];

End[]

EndPackage[]
