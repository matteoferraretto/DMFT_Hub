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
EdMode = "Normal";
U = {-1.0};
Ust = 0.0;
Usec = 0.0
Jph = 0.0;
Jse = 0.0;
\[Mu] = 0;

BathParameters = Flatten[StartingBath["Default", L-1, Norb, EdMode]]
InteractionParameters = Flatten[{U, Ust, Usec, Jph, Jse, \[Mu]}];

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

HImp[L, f, Norb, Sectors, BathParameters, InteractionParameters, EdMode]


(* ::Subtitle:: *)
(*Test the function for Norb = 2*)


ClearAll["Global*`"];
L = 3;
f = 2;
Norb = 2;
EdMode = "Normal";
U = {-1.0, -1.0};
Ust = 0.0;
Usec = 0.0
Jph = 0.0;
Jse = 0.0;
\[Mu] = 0;

BathParameters = Flatten[StartingBath["Default", L-1, Norb, EdMode]]
InteractionParameters = Flatten[{U, Ust, Usec, Jph, Jse, \[Mu]}];

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

HImp[L, f, Norb, Sectors, BathParameters, InteractionParameters, EdMode]



ClearAll["Global*`"];
L = 3;
f = 2;
Norb = 2;
EdMode = "Superc";
U = {-1.0, -1.0};
Ust = 0.0;
Usec = 0.0
Jph = 0.0;
Jse = 0.0;
\[Mu] = 0;

BathParameters = Flatten[StartingBath["Default", L-1, Norb, EdMode]]
InteractionParameters = Flatten[{U, Ust, Usec, Jph, Jse, \[Mu]}];

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

HImp[L, f, Norb, Sectors, BathParameters, InteractionParameters, EdMode]


ClearAll["Global*`"];
L = 2;
f = 2;
Norb = 2;
EdMode = "FullSuperc";
U = {-1.0, -1.0};
Ust = 0.0;
Usec = 0.0
Jph = 0.0;
Jse = 0.0;
\[Mu] = 0;

BathParameters = Flatten@StartingBath["Default", L-1, Norb, "Superc"]

InteractionParameters = Flatten[{U, Ust, Usec, Jph, Jse, \[Mu]}]

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

HImp[L, f, Norb, Sectors, BathParameters, InteractionParameters, EdMode]




