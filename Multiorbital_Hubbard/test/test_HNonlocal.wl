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

HNonlocal[L, f, Norb, Sectors, EdMode];
Dimensions[%]

(* build H *)
H = Dot[Parameters, #]&/@HNonlocal[L, f, Norb, Sectors, EdMode];
(* check that the Hamiltonian has as many blocks as the number of sectors *)
Length[H] == Length[Sectors]
(* check that each block of the Hamiltonian is hermitian *)
HermitianMatrixQ[#]&/@H
(* print one of the blocks *)
H[[6]]//MatrixForm



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
(**)
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

HNonlocal[L, f, Norb, Sectors, EdMode];
Dimensions[%]

(* build H *)
H = Dot[Parameters, #]&/@HNonlocal[L, f, Norb, Sectors, EdMode];
(* check that the Hamiltonian has as many blocks as the number of sectors *)
Length[H] == Length[Sectors]
(* check that each block of the Hamiltonian is hermitian *)
HermitianMatrixQ[#]&/@H
(* print one of the blocks *)
H[[6]]//MatrixForm


(* ::Subtitle:: *)
(*Test performance*)


ClearAll["Global*`"]
L = 5;
f = 2;
Norb = 2;
EdMode = "Normal";

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;
DimSector[L, f, Norb, #, EdMode]&/@QnsSectorList

AbsoluteTiming[
	HNonlocal[L, f, Norb, Sectors, EdMode];
]



