(* ::Package:: *)

ClearAll["Global*`"];
ClearAll["DMFT*`"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*Test the function for Norb = 1*)


ClearAll["Global*`"];
L = 3;
f = 3;
Norb = 1;
EdMode = "Raman";
U = 1.0;
\[Mu] = 0;
h = 2.0 * Table[-1/2.+k/(f-1), {k, 0, f-1}];
Parameters = Flatten[{U, \[Mu], h}]

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

AbsoluteTiming[
	HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode];
]
Dimensions[#]&/@HlocBlocks


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
L = 3;
f = 2;
Norb = 2;
EdMode = "Raman";
U = {1.0, 1.0};
Ust = 0;
JH = 0;
\[Mu] = 0.0;
h = 2.0 * Table[-1/2.+k/(f-1), {k,0,f-1}];
Parameters = Flatten[{U, Ust-2*JH, Ust-3*JH, \[Mu], ConstantArray[h, Norb]}]

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

AbsoluteTiming[
	HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode];
]

Dimensions[Parameters]
Dimensions[HlocBlocks]

Hloc = SparseArray[#]&/@(
	Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@HlocBlocks
	);
	

(*
(Eigenvalues[#]&/@Hloc)[[150]]
(Eigenvalues[#]&/@Hloc)[[90]]
QnsSectorList[[150]]
QnsSectorList[[90]]
IntegerDigits[#,2,L]&/@Sectors[[150]]
IntegerDigits[#,2,L]&/@Sectors[[90]]
*)
Dimensions[Hloc]




