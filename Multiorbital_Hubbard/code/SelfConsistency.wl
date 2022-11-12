(* ::Package:: *)

BeginPackage["SelfConsistency`", {"MyLinearAlgebra`"}]


(* Self consistency *)
WeissField::usage = "WeissField[L_, f_, \[Mu]_, symbols_, z_, EdMode_]"

SelfConsistency::usage = "SelfConsistency[DBethe, \[Mu], Weiss, symbols, StartingParameters, LocalG, zlist, EdMode]"

DMFTError::usage = "DMFTError[Xnew, Xold, EdMode]"


Begin["Private`"]

Print["Package SelfConsistency` loaded successfully."];

(* ANALYTIC EVALUATION OF NONINTERACTING IMPURITY GREEN FUNCTION *)
(* symbolic non-interacting Green function obtained by inverting the full impurity-bath Hamiltonian *)
GreenFunction0[L_, f_,\[Mu]_, symbols_, z_, EdMode_] := Module[
	{e, V, \[CapitalDelta], H, G},
	Which[
		EdMode == "Normal",
		e = symbols[[1;;L-1]];
		V = symbols[[L;;2(L-1)]];
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
		e = symbols[[1;;L-1]];
		V = symbols[[L;;2(L-1)]];
		\[CapitalDelta] = symbols[[2L-1;;3(L-1)]];
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
			InverseElement[SparseArray[z*IdentityMatrix[f*L] - H], {i,j}]
		,{i,1,2},{j,1,2}];(* the IMPURITY part of the Green function is the 2x2 top left block *)
	];
G
];

(* Weiss field: non-interacting inverse of the Green function *)
WeissField[L_, f_, \[Mu]_, symbols_, z_, EdMode_] := With[
	{G0 = GreenFunction0[L, f, \[Mu], symbols, z, EdMode]},
	Which[
		EdMode == "Normal", FullSimplify[1/G0],
		EdMode == "Superc", Inverse[G0]
	]
];

(* Perform minimization according to the self consistency condition *)
SelfConsistency[DBethe_, \[Mu]_, Weiss_, symbols_, z_, IndependentParameters_, LocalG_, \[CapitalSigma]_, zlist_, EdMode_, OptionsPattern[{SelfConsistency, FindMinimum}]] := Module[
	{Lattice = OptionValue[Lattice],
	LFit = Min[Length[zlist], OptionValue[NumberOfFrequencies]],
	weight = OptionValue[FitWeight],
	residue, newparameters, \[Chi], \[Sigma]3 = PauliMatrix[3]},
	(* define the target function to minimize depending on the Lattice and EdMode (if Lattice = Bethe there is a shortcut) *)
	Which[
		EdMode == "Normal" && Lattice == "Bethe",
		\[Chi][symbols] = Mean[Abs[ weight[[;;LFit]] * (
			((Weiss - z - \[Mu])/.{z -> #}&/@Take[zlist, LFit]) + (DBethe^2/4.)*Take[LocalG, LFit]
		)]^2],
	(* ------------------------------------ *)
		EdMode == "Superc" && Lattice == "Bethe",
		\[Chi][symbols] = Mean@First@
			Mean[Abs[
				((Weiss - z - \[Mu]*\[Sigma]3)/.{z -> #}&/@Take[zlist, LFit]) 
				+ (DBethe^2/4.)*Map[Dot[\[Sigma]3, #, \[Sigma]3]&, Take[LocalG, LFit]]
			]^2],
	(* ------------------------------------ *)
		EdMode == "Normal" && Lattice != "Bethe",
		\[Chi][symbols] = Mean[Abs[ weight[[;;LFit]] * (
			(Weiss/.{z -> #}&/@zlist[[;;LFit]]) - (1./LocalG + \[CapitalSigma])[[;;LFit]]
		)]^2],
	(* ------------------------------------ *)
		EdMode == "Superc" && Lattice != "Bethe",
		\[Chi][symbols] = Mean@First@
			Mean[Abs[
				(Weiss/.{z -> #}& /@ zlist[[;;LFit]]) - (Inverse/@(LocalG[[;;LFit]]) + \[CapitalSigma][[;;LFit]])
			]^2]
	];
	Which[
		OptionValue[Minimum] == "Local",
		{residue, newparameters} =
			FindMinimum[
				\[Chi][symbols],
				{symbols, IndependentParameters}\[Transpose],
				Method -> "ConjugateGradient",
				MaxIterations -> OptionValue[MaxIterations],
				AccuracyGoal -> OptionValue[AccuracyGoal]
			];,
	(* -------------------------------------- *)
		OptionValue[Minimum] == "Global",
		{residue, newparameters} =
			NMinimize[
				\[Chi][symbols],
				symbols,
				Method -> OptionValue[Method],
				MaxIterations -> OptionValue[MaxIterations],
				AccuracyGoal -> OptionValue[AccuracyGoal]
			];
	];
	Print["Fit residue: ", residue];
	symbols/.newparameters
];
Options[SelfConsistency] = {Lattice -> "Bethe", LatticeDimension -> 2, NumberOfFrequencies -> 2000, Minimum -> "Local", FitWeight -> ConstantArray[1., 2000]};

(* DMFT error *)
DMFTError[Xnew_, Xold_, EdMode_] := Module[
	{error},
	Which[
		EdMode == "Normal",
		error = Total[Abs[
			Xnew - Xold
		]]/Max[
			Total[Abs[Xnew]], Total[Abs[Xold]]
		],
(* ------------------------------- *)
		EdMode == "Superc",
		error = Mean[{
			Total[Abs[
				Xnew[[All, 1, 1]] - Xold[[All, 1, 1]]
			]]/Max[
				Total[Abs[Xnew[[All, 1, 1]]]], Total[Abs[Xold[[All, 1, 1]]]]		
			],
			Total[Abs[
				Xnew[[All, 1, 2]] - Xold[[All, 1, 2]]
			]]/Max[
				Total[Abs[Xnew[[All, 1, 2]]]], Total[Abs[Xold[[All, 1, 2]]]]		
			]
		}]
	];
	error
];

End[]

EndPackage[]
