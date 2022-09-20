(* ::Package:: *)

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


L = 1;
f = 2;
Norb = 2;
EdMode = "Normal";
U1 = -1.0;
U2 = -1.0;
Ust = -1.0;
JH = 0.1;
Parameters = {U1, U2, Ust-2*JH, Ust-3*JH, 0};

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

AbsoluteTiming[
	Hlocal = Dot[Parameters, #]&/@ImpHLocalNormal[L,f,Norb,Sectors];
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
