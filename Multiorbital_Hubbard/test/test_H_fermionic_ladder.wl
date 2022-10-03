(* ::Package:: *)

ClearAll["Global*`"];
ClearAll["DMFT*`"];

FolderPath = NotebookDirectory[];
<<(FolderPath<>"/../Multiorbital_Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*Test the function for Norb = 1*)


ClearAll["Global*`"];
L = 2;
f = 2;
Norb = 1;
EdMode = "Raman";

e = 0.2;(* electric field *)
t = -1.0;(* real hopping *)
\[CapitalOmega] = 0.5;(* Raman hopping *)
\[Gamma] = Pi/2;(* effective flux *)
\[Phi] = Pi/10;(* phase of real hopping *)


V = -e*Table[j - (L+1)/2., {j, L}](* potential *)
Parameters = Join[
	Flatten@ConstantArray[V, f],
	ConstantArray[t, f*(L-1)],
	ConstantArray[\[CapitalOmega], (f-1)*L]
](* OBCs in both directions *)
Length[Parameters]

QnsSectorList = SectorList[L, f, Norb, EdMode];
Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;

HalfFilledSector = BuildSector[L, f, Norb, {L}, EdMode];
IntegerDigits[#,2,L]&/@HalfFilledSector;

Hblocks = First @ Hnonint[L, f, Norb, {HalfFilledSector}, EdMode, RealPBC -> False, SyntheticPBC -> False, RealPhase -> ConstantArray[\[Phi], f], SyntheticPhase -> \[Gamma]];
Length[Hblocks]

H = Sum[
		Parameters[[i]]*Hblocks[[i]], 
	{i, Length@Parameters}]
	
HermitianMatrixQ[H]
H//MatrixForm
		
Sort@Eigenvalues[H]


(* ::Subtitle:: *)
(**)
(*Performance test*)


ClearAll["Global*`"];
L = 6;
f = 3;
Norb = 1;
EdMode = "Raman";

e = 0.2;(* electric field *)
t = -1.0;(* real hopping *)
\[CapitalOmega] = 0.5;(* Raman hopping *)
\[Gamma] = Pi/2;(* effective flux *)
\[Phi] = Pi/10;(* phase of real hopping *)

V = -e*Table[j - (L+1)/2., {j, L}];(* potential *)
Parameters = Join[
	Flatten@ConstantArray[V, f],
	ConstantArray[t, f*(L-1)],
	ConstantArray[\[CapitalOmega], (f-1)*L]
];(* OBCs in both directions *)


AbsoluteTiming[
	QnsSectorList = SectorList[L, f, Norb, EdMode];
	Sectors = BuildSector[L, f, Norb, #, EdMode]&/@QnsSectorList;
	HalfFilledSector = BuildSector[L, f, Norb, {L}, EdMode];
	Print["half filled sector dimension: ", Length@HalfFilledSector];
	Hblocks = First @ Hnonint[L, f, Norb, {HalfFilledSector}, EdMode, RealPBC -> False, SyntheticPBC -> False, RealPhase -> ConstantArray[\[Phi], f], SyntheticPhase -> \[Gamma]];
	H = Sum[
		Parameters[[i]]*Hblocks[[i]], 
	{i, Length@Parameters}];
]



