(* ::Package:: *)

ClearAll["Global`*"];
ClearAll["DMFT`*"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*Check that it works with Norb = 1*)


ClearAll["Global`*"]
L = 3;
f = 2;
Norb = 1;
EdMode = "Raman";

Parameters = StartingBath[L, f, Norb, "Default", EdMode]
Parameters = Flatten[Parameters]
Length[Parameters]

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList
Length[QnsSectorList]

HnonlocBlocks = HNonlocalRaman[L, f, Norb, Sectors, EdMode];

Table[
	Dimensions[HnonlocBlocks[[i]]],
	{i, Length@QnsSectorList}]

Hnonloc = SparseArray[#]&/@(
	Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@HnonlocBlocks
	)
Dimensions[Hnonloc];
(* build H 
H = Dot[Parameters, #]&/@HNonlocal[L, f, Norb, Sectors, EdMode]; *)
(* check that the Hamiltonian has as many blocks as the number of sectors *)
Dimensions[HnonlocBlocks]
ByteCount[Hnonloc]
(* check that each block of the Hamiltonian is hermitian *)
HermitianMatrixQ[#]&/@Hnonloc



(* ::Subtitle:: *)
(*Test performance with large matrices*)


ClearAll["Global`*"]
L = 8;
f = 2;
Norb = 1;
EdMode = "Raman";

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;
Length[#]&/@Sectors

AbsoluteTiming[
	HNonlocalRaman[L, f, Norb, Sectors, EdMode];
]



(* ::Subtitle:: *)
(*Check that it works with Norb = 2*)


ClearAll["Global`*"]
L = 3;
f = 2;
Norb = 2;
EdMode = "Raman";

Parameters = StartingBath[L, f, Norb, "Default", EdMode];
Parameters = Flatten[Parameters]

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

(* build H *)

HnonlocBlocks = HNonlocalRaman[L, f, Norb, Sectors, EdMode];
Hnonloc = SparseArray[#]&/@(
	Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@HnonlocBlocks
	)

(* check dimensions *)
Dimensions[HnonlocBlocks]
ByteCount[Hnonloc]
(* check that the Hamiltonian has as many blocks as the number of sectors *)
Length[Hnonloc] == Length[Sectors]
(* check that each block of the Hamiltonian is hermitian *)
HermitianMatrixQ[#]&/@Hnonloc



(* ::Subtitle:: *)
(*Test performance*)


ClearAll["Global`*"]
L = 5;
f = 2;
Norb = 2;
EdMode = "Raman";

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;
Length[#]&/@Sectors

AbsoluteTiming[
	HNonlocalRaman[L, f, Norb, Sectors, EdMode];
]



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
				Print["flag=",flag,". orb=", orb,". {\[Rho],\[Sigma]}=",{\[Rho],\[Sigma]},". j=",j];
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
				Print["flag=",flag,". orb=", orb,". {\[Rho],\[Sigma]}=",{\[Rho],\[Sigma]},". j=",j];
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
		,{flag, {"Bath", "Hopping"}}, {orb, 1, Norb}, {j, OptionValue[Nimp]+1, L}, {\[Rho], 1, f}, {\[Sigma], 1, f}];
		(*AppendTo[H, Hsector];*)
	,{\[Psi], Sectors}];
	H
];
Options[HNonlocalRaman] = {Nimp -> 1};


HNonlocalRaman[L, f, Norb, Sectors, "Raman"]



