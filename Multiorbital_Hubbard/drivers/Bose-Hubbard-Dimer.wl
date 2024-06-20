(* ::Package:: *)

(* generates the list of bosonic Fock states with exactly nbosons on a lattice with L sites *)
BosonicBASIS[L_, nbosons_]:= SortBy[
	Flatten[
		Permutations[#] &/@ (
			PadRight[#, L] &/@ IntegerPartitions[nbosons, L]
		),1]
	, First];

(* gives True if bosonic hopping from site j to site i is possible on a given state. *)
BosonicHopQ = Compile[{
	{i,_Integer}, {j,_Integer}, {state,_Integer,1}
	},
	state[[j]]!=0, 
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True, RuntimeOptions->"Speed"
];

(* select states on the Hilber space for which hopping j->i is possible *)
BosonicHopSelect[i_, j_, stateList_] := With[
	{criteria = BosonicHopQ[i,j,stateList]},
	Pick[stateList, criteria]
];

(* actually perform a hopping j->i on a given state *)
BosonicHop = Compile[{
	{i, _Integer}, {j, _Integer}, {state, _Integer,1}
	},
	Module[{newstate = state},
		newstate[[i]] = state[[i]] + 1;
		newstate[[j]] = state[[j]] - 1;
		newstate
	],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True, RuntimeOptions->"Speed"
];

(* return the coefficient of a bosonc hopping j->i on a given state *)
BosonicHopCoefficient = Compile[{
	{i,_Integer}, {j,_Integer}, {state,_Integer,1}
	},
	Sqrt[1.0*(state[[i]]+1)*state[[j]]],
	RuntimeAttributes->{Listable}, CompilationTarget->"C"
];

(* count the number of bosons in site j for a given state *)
BosonicN = Compile[{
	{i,_Integer}, {state,_Integer,1}},
	state[[i]],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True, RuntimeOptions->"Speed"
];

(* neighbor of site j on a L-sites ring (already in library!) *)
Neighbor = Compile[{
	{L,_Integer}, {j,_Integer}
	},
	Mod[j,L]+1,
	CompilationTarget->"C"
];

(* non-interacting blocks of the Bosonic Hamiltonian on a 1d lattice *)
BosonicHnonint[L_, Sectors_, OptionsPattern[]] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,dimeff, rules,dispatch,cols,rows,pos,num, coeff},
	H = Last @ Last @ Reap[
		Do[
			Hsector = {};
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			Hsector = Last @ Last @ Reap[
				Do[
					Hblock = SparseArray[{}, {dim,dim}];
					Which[
						flag == "Potential",
						num = BosonicN[j, \[Psi]];(*local density*)
						Hblock += SparseArray @ DiagonalMatrix[num];
						Sow[Hblock, "Hblock"];,
					(* --------------------------------------------------------------- *)
						flag == "Hopping",
						If[j == L && !OptionValue["RealPBC"], Continue[];];
						\[Psi]1 = BosonicHopSelect[j, Neighbor[L,j], \[Psi]];
						dimeff = Length[\[Psi]1];
						If[dimeff == 0, Continue[];];
						\[Chi] = BosonicHop[j, Neighbor[L,j], \[Psi]1];
						coeff = BosonicHopCoefficient[j, Neighbor[L,j], \[Psi]1];
						rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						Hblock += SparseArray[pos -> coeff, {dim,dim}];
						Hblock = Hblock + Hblock\[ConjugateTranspose];
						Sow[Hblock, "Hblock"];
					];
				, {flag, {"Potential", "Hopping"}}, {j,1,L}];
			, "Hblock"];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"];
	H
];
Options[BosonicHnonint] = {"RealPBC" -> False, "RealPhase" -> 0.0};

(* matrix associated to hopping j -> i written in the Fock basis *)
CoherenceVisibilityMatrix[L_, Sectors_, i_, j_] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,dimeff, rules,dispatch,cols,rows,pos,coeff,num},
	If[i>L || j>L || i<1 || j<1, Print["Error, i and j should be within the range [1, L]."]];
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
			dispatch = Dispatch[rules];
			\[Psi]1 = BosonicHopSelect[i, j, \[Psi]];
			Hsector = SparseArray[{}, {dim, dim}];
			dimeff = Length[\[Psi]1];
			If[dimeff == 0, Continue[];];
			\[Chi] = BosonicHop[i, j, \[Psi]1];
			coeff = BosonicHopCoefficient[i, j, \[Psi]1];
			rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
			Hsector = SparseArray[pos -> coeff, {dim,dim}];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"];
	H
];

(* interacting Hamiltonian *)
BosonicHint[L_, Sectors_] := Module[
	{H, num, Hsector, Hblock, dim},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			Hsector = SparseArray[{},{dim,dim}];
			Do[
				num = BosonicN[j, \[Psi]];
				Hsector += SparseArray @ DiagonalMatrix[0.5*num*(num-1)];
			, {j, 1, L}];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"]
];

DensityDifferenceToPowerK[power_, L_, Sectors_] := Module[
	{H, num, Hsector, Hblock, dim},
	H = Last @ Last @ Reap[
		Do[
			dim = Length[\[Psi]];
			Hsector = SparseArray[{}, {dim, dim}];
			num = (BosonicN[1, \[Psi]] - BosonicN[2,\[Psi]])^power;
			Hsector += SparseArray @ DiagonalMatrix[num];
			Sow[Hsector, "Hsector"];
		, {\[Psi], Sectors}];
	, "Hsector"]
];

BoseHubbardHamiltonian[L_, Sectors_, \[Epsilon]list_,J_, U_] := Module[
	{H, parameters},
	parameters = Join[\[Epsilon]list, ConstantArray[-J, L-1]];
	H = Sum[
		parameters[[n]] * #[[n]]
	, {n, Length[parameters]}] &/@ BosonicHnonint[L, Sectors] + U*BosonicHint[L, Sectors]
];

Eigs[H_] := SortBy[
	Eigensystem[H]\[Transpose]
, First]\[Transpose];

(* returns the density matrix written in the Fock basis *)
DensityMatrix[eigs_, T_] := Module[
	{dim, evals = eigs[[1]], evecs = eigs[[2]], d, p, Z, threshold = 10^(-16)},
	(* Hilbert space dimension *)
	dim = Length[evals];
	(* select values of \[Beta](E-E_0) that are not too large, and exponentiate it  *)
	evals = Select[(evals - evals[[1]])/T, #<-Log[threshold]&];
	evals = Exp[-evals];
	(* compute the partition function *)
	Z = Total[evals];
	(* density matrix in the basis of the eigenstates *)
	d = DiagonalMatrix[
		PadRight[evals/Z, dim]
	];
	(* rotation matrix *)
	p = evecs\[Transpose];
	p . d . Inverse[p]
];

EntanglementEntropy = Compile[{
	{\[Rho], _Complex, 2}
	},
	With[
		{p = Table[\[Rho][[j,j]], {j, Length[\[Rho]]}]},
		-p . Log[2,p]
	],
	RuntimeAttributes->{Listable}, CompilationTarget->"C", Parallelization->True
];

FisherInformation[\[Rho]_, nbosons_] := Module[
	{L = 2, sector, densdiffsquare, densdiff},
	sector = BosonicBASIS[L, nbosons];
	densdiffsquare = DensityDifferenceToPowerK[2, L, {sector}][[1]];
	densdiff = DensityDifferenceToPowerK[1, L, {sector}][[1]];
	(Tr[\[Rho] . densdiffsquare] - Tr[\[Rho] . densdiff]^2)/(nbosons^2)
];
