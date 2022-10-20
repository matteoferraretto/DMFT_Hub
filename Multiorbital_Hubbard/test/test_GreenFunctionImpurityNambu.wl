(* ::Package:: *)

ClearAll["Global`*"];
ClearAll["DMFT`*"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*Test for the non interacting chain*)


ClearAll["Global`*"];
L = 4;
f = 2;
Norb = 1;
EdMode = "Superc";
DegeneracyThreshold = 10^(-9);
\[Mu] = 0;
Parameters = StartingBath[L, f, Norb, "Default", EdMode]
FlatParameters = Flatten[Parameters]

(* reshape parameters to a list suitable for computation of non interacting G.F. *)
ReshapeParameters[L_, f_, \[Sigma]_, orb_, Parameters_, EdMode_] := Which[
	EdMode == "Normal",
	Thread[
		Join[
			Table[Symbol["e"<>ToString[i]], {i, L-1}],
			Table[Symbol["V"<>ToString[i]], {i, L-1}]
		] -> 
		Join[
			Parameters[[1]][[f*(orb-1)+\[Sigma]]],(* e_1,\[Sigma],orb , e_2,\[Sigma],orb , ... *)
			Parameters[[2]][[f*(orb-1)+\[Sigma]]](* V_1,\[Sigma],orb , V_2,\[Sigma],orb , ... *)
		]
	],
	(* --------------------------- *)
	EdMode == "Superc",
	Thread[
		Join[
			Table[Symbol["e"<>ToString[i]], {i, L-1}],
			Table[Symbol["V"<>ToString[i]], {i, L-1}],
			Table[Symbol["\[CapitalDelta]"<>ToString[i]], {i, L-1}]
		] ->
		Join[
			Parameters[[1]][[f*(orb-1)+\[Sigma]]],(* e_1,\[Sigma],orb , e_2,\[Sigma],orb , ... *)
			Parameters[[2]][[f*(orb-1)+\[Sigma]]],(* V_1,\[Sigma],orb , V_2,\[Sigma],orb , ... *)
			Parameters[[3]][[orb]](* \[CapitalDelta]_1,\[Sigma],orb , \[CapitalDelta]_2,\[Sigma],orb , ... *)
		]
	]
];



QnsSectorList = SectorList[L, f, Norb, EdMode]
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;
SectorsDispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,QnsSectorList],1]];

Hblocks = HNonlocal[L, f, Norb, Sectors, EdMode];
Hsectors = SparseArray[#]&/@(
	Sum[FlatParameters[[i]]*#[[i]],{i,1,Length@FlatParameters}]&/@Hblocks
);

(*energies = Eigenvalues[#]&/@Hsectors
eigenstates = Eigenvectors[#]&/@H;*)
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


InverseElement[m_, {i_,j_}] := (-1)^(i+j)Det[Drop[m,{j},{i}]]/Det[m];

GreenFunction0[L_,f_,EdMode_]:=Module[{
e=Table[Symbol["e"<>ToString[i]],{i,L-1}],
V=Table[Symbol["V"<>ToString[i]],{i,L-1}],
\[CapitalDelta]=Table[Symbol["\[CapitalDelta]"<>ToString[i]],{i,L-1}],
H, G},
Which[
EdMode == "Normal",
(* the spinor is (d, c_1, c_2, ...) *)
H = SparseArray[{
{i_,j_}/;(i==1&&j>1):>V[[j-1]]
},
{L,L}];
H = H + H\[Transpose];
H += SparseArray[{
{1,1}->-\[Mu],
{i_,i_}/;(i>1):>e[[i-1]]
},
{L,L}];
G = InverseElement[SparseArray[z*IdentityMatrix[L]-H], {1,1}];,
(* ------------------------ *)
EdMode == "Superc",
(* the spinor is (d_up, ddg_dw, c_1up, cdg_1dw, c_2up, cdg_2dw, ... *)
H = SparseArray[{
{i_,j_}/;(j==i+1&&i>2&&Mod[i,f]==1):>\[CapitalDelta][[Quotient[i,f]]],
{i_,j_}/;(i==1&&j>2&&Mod[j,f]==1):>V[[Quotient[j,f]]],
{i_,j_}/;(i==2&&j>2&&Mod[j,f]==0):>-V[[Quotient[j-1,f]]]
},
{f*L, f*L}];
H = H + H\[Transpose];
H += SparseArray[{
{1, 1}->-\[Mu], {2,2}->\[Mu],
{i_,i_}/;(i>2&&Mod[i,f]==1):>e[[Quotient[i-1,f]]],
{i_,i_}/;(i>2&&Mod[i,f]==0):>-e[[Quotient[(i-2),f]]]
},
{f*L,f*L}
];
G = Table[
InverseElement[SparseArray[z*IdentityMatrix[f*L]-H], {i,j}]
,{i,1,2},{j,1,2}];(* the IMPURITY part of the Green function is the 2x2 top left block *)
];
G
];

GF0[z] = GreenFunction0[L,f,EdMode]/.ReshapeParameters[L, f, 1, 1, Parameters, EdMode]


d\[Omega] = 0.005; \[Eta] = 0.1;
zlist = Table[n*d\[Omega], {n, -800,800}] + \[Eta]*I;

AbsoluteTiming[
	gfLanczos = GreenFunctionImpurityNambu[L, f, Norb, 1, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist];
]
AbsoluteTiming[
	gfAnalytic = (GF0[z]/.{z->#})&/@zlist;
]

ListPlot[{
	{Re[zlist], -(1/Pi)*Im[gfLanczos[[All,1,1]]]}\[Transpose],
	{Re[zlist], -(1/Pi)*Im[gfAnalytic[[All,1,1]]]}\[Transpose]
	},
	PlotStyle->{Thick, Dashing[.05]},
	PlotRange->All,
	Joined->True,
	AxesLabel->{"\[Omega]","DoS"}]
	
ListPlot[{
	{Re[zlist], -(1/Pi)*Im[gfLanczos[[All,1,2]]]}\[Transpose],
	{Re[zlist], -(1/Pi)*Im[gfAnalytic[[All,1,2]]]}\[Transpose]
	},
	PlotStyle->{Thick, Dashing[.05]},
	PlotRange->All,
	Joined->True,
	AxesLabel->{"\[Omega]","DoS"}]
	
ListPlot[{
	{Re[zlist], -(1/Pi)*Im[gfLanczos[[All,2,2]]]}\[Transpose],
	{Re[zlist], -(1/Pi)*Im[gfAnalytic[[All,2,2]]]}\[Transpose]
	},
	PlotStyle->{Thick, Dashing[.05]},
	PlotRange->All,
	Joined->True,
	AxesLabel->{"\[Omega]","DoS"}]
	
ListPlot[{
	{Re[zlist], -(1/Pi)*Im[gfLanczos[[All,2,1]]]}\[Transpose],
	{Re[zlist], -(1/Pi)*Im[gfAnalytic[[All,2,1]]]}\[Transpose]
	},
	PlotStyle->{Thick, Dashing[.05]},
	PlotRange->All,
	Joined->True,
	AxesLabel->{"\[Omega]","DoS"}]
	
(* total area under the spectral function. It should be around 1, not exactly 1 due to the presence of broadening \[Eta] *)
Total[
	-(1/Pi)*Im[gfLanczos]*d\[Omega]
]









