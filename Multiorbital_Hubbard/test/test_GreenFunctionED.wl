(* ::Package:: *)

ClearAll["Global`*"];
ClearAll["DMFT`*"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*Test the function for Norb = 1*)


ClearAll["Global`*"];
L = 4;
f = 1;
Norb = 1;
EdMode = "Normal";

e = 0;(* electric field *)
t = 1.0;(* real hopping *)
nParticles = 2;(* number of particles *)
T = 0.1;(* temperature *)

realPBC = True;


(* local potential *)
V = -e*Table[j - (L+1)/2., {j, L}];
(* List of physical parameters in correct order *)
Parameters = Join[
	Flatten@ConstantArray[V, f],
	ConstantArray[-t, If[realPBC, f*L, (*else*)f*(L-1)]]
]
Length[Parameters]

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList
SectorsDispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,QnsSectorList],1]]

Hblocks = Hnonint[L, f, Norb, Sectors, EdMode, RealPBC -> realPBC, RealPhase -> {0}];
Dimensions[Hblocks]

H = SparseArray[#]&/@(
	Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@Hblocks
)

(*energies = Eigenvalues[#]&/@H;
eigenstates = Eigenvectors[#]&/@H;*)
eigs = Eigensystem[#]&/@H;
energies = eigs[[All,1]];
eigenstates = eigs[[All,2]];


d\[Omega] = 0.005; \[Eta] = 0.1;
zlist = Table[n*d\[Omega], {n, -800,800}] + \[Eta]*I;
\[Epsilon][k_] := -2*t*Cos[k];
kvalues = Table[2Pi*n/L, {n, 0, L-1}]

gf = GreenFunctionED[L, f, Norb, {1,1}, 1, 1, Sectors, QnsSectorList, eigs, 0.01, zlist, EdMode, NormalizedFunction->True];

ListPlot[
	{Re[zlist], -(1/Pi)*Im[gf]}\[Transpose],
	PlotRange->All,
	AxesLabel->{"\[Omega]","DoS"}]
	
(* ----------------------------------------- *)
Show[
	Plot[\[Epsilon][k], {k,0,2Pi},AxesLabel->{"k","\[Omega]"}],
	ListPlot[{kvalues, \[Epsilon]/@kvalues}\[Transpose],PlotStyle->{Blue,PointSize[.02]}]
]


(* total area under the spectral function. It should be around 1, not exactly 1 due to the presence of broadening \[Eta] *)
Total[
	-(1/Pi)*Im[gf]*d\[Omega]
]
