(* ::Package:: *)

ClearAll["Global`*"];
ClearAll["DMFT`*"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*Test for the non interacting chain*)


ClearAll["Global`*"];
L = 6;
f = 1;
Norb = 1;
EdMode = "Normal";
DegeneracyThreshold = 10^(-9);

e = 0;(* electric field *)
t = 1.0;(* real hopping *)
T = 0.1;(* temperature *)

realPBC = True;


(* local potential *)
V = -e*Table[j - (L+1)/2., {j, L}];
(* List of physical parameters in correct order *)
Parameters = Join[
	Flatten@ConstantArray[V, f],
	ConstantArray[-t, If[realPBC, f*L, (*else*)f*(L-1)]]
];

QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;
SectorsDispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,QnsSectorList],1]];

Hblocks = Hnonint[L, f, Norb, Sectors, EdMode, RealPBC -> realPBC, RealPhase -> {0}];
Hsectors = SparseArray[#]&/@(
	Sum[Parameters[[i]]*#[[i]],{i,1,Length@Parameters}]&/@Hblocks
);

(*energies = Eigenvalues[#]&/@Hsectors
eigenstates = Eigenvectors[#]&/@H;*)
eigs = Eigensystem[#]&/@Hsectors;
{EgsSectorList, GsSectorList} = Map[Sort[Transpose[Eigensystem[#]]][[1]]&,Hsectors]\[Transpose];(*find the ground state for each sector*)
EgsSectorList = Flatten@EgsSectorList(*correctly reshape the list *)
GsSectorList = Replace[GsSectorList,{x_List}:>x,{0,-2}](*correctly reshape the list*)

(* Compute the ground state *)
Egs = Min@EgsSectorList(*ground state energy (lowest of all the sectors)*)
GsSectorIndex =
	Flatten@Position[EgsSectorList,
			_?((Abs[#-Egs] < DegeneracyThreshold)&)
		];(*sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy*)
DegeneracyWarning = If[Length[GsSectorIndex]>1, True, (*else*) False];(*is True if the ground state is degenerate, False otherwise*)
GsQns = QnsSectorList[[GsSectorIndex]][[1]]

Gs = GsSectorList[[GsQns/.SectorsDispatch]]


d\[Omega] = 0.005; \[Eta] = 0.1;
zlist = Table[n*d\[Omega], {n, -800,800}] + \[Eta]*I;
\[Epsilon][k_] := -2*t*Cos[k];
kvalues = Table[2Pi*n/L, {n, 0, L-1}]

AbsoluteTiming[
	gfLanczos = GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist];
]

AbsoluteTiming[
	gfED = GreenFunctionED[L, f, Norb, {1,1}, 1, 1, Sectors, QnsSectorList, eigs, T, zlist, EdMode, NormalizedFunction->True];
]

ListPlot[{
	{Re[zlist], -(1/Pi)*Im[gfLanczos]}\[Transpose],
	{Re[zlist], -(1/Pi)*Im[gfED]}\[Transpose]
	},
	PlotStyle->{Thick, Dashing[.02]},
	PlotRange->All,
	Joined->True,
	AxesLabel->{"\[Omega]","DoS"}]
	
(* ----------------------------------------- *)
Show[
	Plot[\[Epsilon][k], {k,0,2Pi},AxesLabel->{"k","\[Omega]"}],
	ListPlot[{kvalues, \[Epsilon]/@kvalues}\[Transpose],PlotStyle->{Blue,PointSize[.02]}]
]

(* total area under the spectral function. It should be around 1, not exactly 1 due to the presence of broadening \[Eta] *)
Total[
	-(1/Pi)*Im[gfLanczos]*d\[Omega]
]
Total[
	-(1/Pi)*Im[gfED]*d\[Omega]
]



