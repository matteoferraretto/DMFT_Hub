(* ::Package:: *)

ClearAll["Global*`"];
ClearAll["DMFT*`"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*Check that it works with Norb = 1*)


ClearAll["Global*`"]
L = 3;
f = 2;
Norb = 1;
EdMode = "Superc";

Parameters = StartingBath["Default", L-1, Norb, EdMode];
Parameters = Flatten[Parameters]

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

HnonlocBlocks = HNonlocal[L, f, Norb, Sectors, EdMode];

Hnonloc = SparseArray[#]&/@(
	Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@HnonlocBlocks
	)

(* build H 
H = Dot[Parameters, #]&/@HNonlocal[L, f, Norb, Sectors, EdMode]; *)
(* check that the Hamiltonian has as many blocks as the number of sectors *)
Dimensions[HnonlocBlocks]
ByteCount[Hnonloc]
(* check that each block of the Hamiltonian is hermitian *)
HermitianMatrixQ[#]&/@Hnonloc



(* ::Subtitle:: *)
(*Test performance with large matrices*)


ClearAll["Global*`"]
L = 8;
f = 2;
Norb = 1;
EdMode = "Superc";

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;
DimSector[L, f, Norb, #, EdMode]&/@QnsSectorList

AbsoluteTiming[
	HNonlocal[L, f, Norb, Sectors, EdMode];
]


(* ::Subtitle:: *)
(*Check that it works with Norb = 2*)


ClearAll["Global*`"]
L = 3;
f = 2;
Norb = 2;
EdMode = "Normal";

Parameters = StartingBath["Default", L-1, Norb, EdMode];
Parameters = Flatten[Parameters]

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

(* build H *)
HnonlocBlocks = HNonlocal[L, f, Norb, Sectors, EdMode];
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


ClearAll["Global*`"]
L = 2;
f = 2;
Norb = 2;
EdMode = "FullSuperc";

Parameters = StartingBath["Default", L-1, Norb, EdMode]
Parameters = Flatten[Parameters];

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

(* build H *)
HnonlocBlocks = HNonlocal[L, f, Norb, Sectors, EdMode];
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


ClearAll["Global*`"]
L = 4;
f = 2;
Norb = 2;
EdMode = "FullSuperc";

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;
DimSector[L, f, Norb, #, EdMode]&/@QnsSectorList

AbsoluteTiming[
	HNonlocal[L, f, Norb, Sectors, EdMode];
]




