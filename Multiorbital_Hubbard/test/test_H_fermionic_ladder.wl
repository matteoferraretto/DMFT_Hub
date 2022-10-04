(* ::Package:: *)

ClearAll["Global*`"];
ClearAll["DMFT*`"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*Test the function for Norb = 1*)


ClearAll["Global`*"];
L = 2;
f = 2;
Norb = 1;
EdMode = "Raman";

e = 0.2;(* electric field *)
t = -1.0;(* real hopping *)
\[CapitalOmega] = 0.5;(* Raman hopping *)
\[Gamma] = Pi/2;(* effective flux *)
\[Phi] = Pi/10;(* phase of real hopping *)
U = 2.0;(* Hubbard interaction *)
\[Mu] = 0;(* chemical potential term *)

realPBC = True;
syntheticPBC = True;


(* local potential *)
V = -e*Table[j - (L+1)/2., {j, L}]
(* List of physical parameters in correct order *)
Parameters = Join[
	Flatten@ConstantArray[V, f],
	ConstantArray[t, If[realPBC, f*L, (*else*)f*(L-1)]],
	ConstantArray[\[CapitalOmega], If[syntheticPBC, f*L, (*else*)(f-1)*L]],
	ConstantArray[U, Norb],
	{\[Mu]}
]
Length[Parameters]

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

HalfFilledSector = BuildSector[L, f, Norb, {L}, EdMode];
IntegerDigits[#,2,L]&/@HalfFilledSector

Hblocks = Join[
	First @ Hnonint[L, f, Norb, {HalfFilledSector}, EdMode, RealPBC -> realPBC, SyntheticPBC -> syntheticPBC, RealPhase -> ConstantArray[\[Phi], f], SyntheticPhase -> \[Gamma]],
	First @ HLocal[L, f, Norb, {HalfFilledSector}, EdMode, Nimp -> L]
];
Length[Hblocks]

H = Sum[
		Parameters[[i]]*Hblocks[[i]], 
	{i, Length@Parameters}]
	
HermitianMatrixQ[H]
H//MatrixForm
		
Sort@Eigenvalues[H]


(* ::Subtitle:: *)
(*Performance test*)


ClearAll["Global`*"];
L = 6;
f = 3;
Norb = 1;
EdMode = "Raman";

e = 0.2;(* electric field *)
t = -1.0;(* real hopping *)
\[CapitalOmega] = 0.5;(* Raman hopping *)
\[Gamma] = Pi/2;(* effective flux *)
\[Phi] = Pi/10;(* phase of real hopping *)
U = 2.0;(* Hubbard interaction *)
\[Mu] = 0;(* chemical potential term *)
realPBC = True;
syntheticPBC = True;

(* local potential *)
V = -e*Table[j - (L+1)/2., {j, L}];
(* List of physical parameters in correct order *)
Parameters = Join[
	Flatten@ConstantArray[V, f],
	ConstantArray[t, If[realPBC, f*L, (*else*)f*(L-1)]],
	ConstantArray[\[CapitalOmega], If[syntheticPBC, f*L, (*else*)(f-1)*L]],
	ConstantArray[U, Norb],
	{\[Mu]}
];

AbsoluteTiming[
	QnsSectorList = SectorList[L, f, Norb, EdMode];
	Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;
	HalfFilledSector = BuildSector[L, f, Norb, {L}, EdMode];
	Print["half filled sector dimension: ", Length@HalfFilledSector];
	Hblocks = Join[
		First @ Hnonint[L, f, Norb, {HalfFilledSector}, EdMode, RealPBC -> realPBC, SyntheticPBC -> syntheticPBC, RealPhase -> ConstantArray[\[Phi], f], SyntheticPhase -> \[Gamma]],
		First @ HLocal[L, f, Norb, {HalfFilledSector}, EdMode, Nimp -> L]
	];
	H = Sum[
		Parameters[[i]]*Hblocks[[i]], 
	{i, Length@Parameters}];
]



