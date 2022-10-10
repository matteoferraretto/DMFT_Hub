(* ::Package:: *)

ClearAll["Global*`"];
ClearAll["DMFT*`"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*Test the function for Norb = 1*)


ClearAll["Global*`"];
L = 6;
f = 2;
Norb = 1;
EdMode = "Superc";
U = -1.0;
\[Mu] = 0;
Parameters = {U, \[Mu]};

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

AbsoluteTiming[
	HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode];
	Dimensions[%]
]

Hloc = SparseArray[#]&/@(
	Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@HlocBlocks
	)

Dimensions[Hloc]
ByteCount[Hloc]

(*
AbsoluteTiming[
	Hlocal = Dot[Parameters, #]&/@ImpHLocalDecoupled[L,f,Norb,Sectors];
]

Print["list of energies:"];
EgsList = Eigenvalues[#]&/@Hlocal
Print["min energy:"];
Egs = Min@Flatten@EgsList

Print["index of sector where the ground state is located:"];
gspos = Flatten@Position[
	Flatten[
		Eigenvalues[#, 1]&/@Hlocal
	], Egs]

Part[QnsSectorList, #]&/@gspos
*)



(* ::Subtitle:: *)
(*Test the function for Norb = 2*)


ClearAll["Global*`"];
L = 2;
f = 2;
Norb = 2;
EdMode = "Normal";
U1 = -1.0;
U2 = -1.0;
Ust = -1.0;
JH = 0.1;
\[Mu] = 0.0;
Parameters = {U1, U2, Ust-2*JH, Ust-3*JH, \[Mu]};

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

AbsoluteTiming[
	HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode];
	Dimensions[%]
]

Hloc = SparseArray[#]&/@(
	Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@HlocBlocks
	)

Dimensions[Hloc]
ByteCount[Hloc]



ClearAll["Global*`"];
L = 2;
f = 2;
Norb = 2;
EdMode = "Superc";
U1 = -1.0;
U2 = -1.0;
Ust = -1.0;
JH = 0.1;
Jph = 0.2;
\[Mu] = 0.0;
Parameters = {U1, U2, Ust-2*JH, Ust-3*JH, Jph, \[Mu]};

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

AbsoluteTiming[
	HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode];
	Dimensions[%]
]

Hloc = SparseArray[#]&/@(
		Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@HlocBlocks
	)

Dimensions[Hloc]
ByteCount[Hloc]


ClearAll["Global*`"];
L = 2;
f = 2;
Norb = 2;
EdMode = "InterorbNormal";
U1 = -1.0;
U2 = -1.0;
Ust = -1.0;
JH = 0.1;
Jph = 0.2;
Jse = 0.0;
\[Mu] = 0.0;
Parameters = {U1, U2, Ust-2*JH, Ust-3*JH, Jph, Jse, \[Mu]};
Length@Parameters

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

AbsoluteTiming[
	HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode];
	Dimensions[HlocBlocks]
]

Hloc = SparseArray[#]&/@(
		Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@HlocBlocks
	)

Dimensions[Hloc]
ByteCount[Hloc]




