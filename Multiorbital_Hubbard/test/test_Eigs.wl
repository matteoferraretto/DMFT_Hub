(* ::Package:: *)

Eigs[H_, OptionsPattern[]] := Module[
	{values, vectors, eigs, 
	dim = Length[H], 
	T = OptionValue["Temperature"], 
	\[Epsilon]deg = OptionValue["DegeneracyThreshold"], 
	\[Epsilon]temp = - OptionValue["Temperature"] * Log[OptionValue["BoltzmannThreshold"]]
	},
	If[dim <= OptionValue["MinLanczosDim"],
		(* if the matrix is small, just compute the full spectrum *)
		eigs = Eigensystem[H];
		eigs = SortBy[eigs\[Transpose], First]\[Transpose];
		values = eigs[[1]];
		(* return a suitable number of states depending on temperature *)
		If[T == 0,
			(* return just the ground state, taking into account possible degeneracies *)
			Return[
				Take[#,
					Count[
						values, x_/;(Abs[x - values[[1]]] < \[Epsilon]deg)
					]
				]&/@eigs
			],
		(* else, if T != 0 *)
			Return[
				Take[#,
					Count[
						values, x_/;((Abs[x - values[[1]]] < \[Epsilon]temp) || (Abs[x - values[[1]]] < \[Epsilon]deg))
					]
				]&/@eigs
			]
		],
	(* else *)
		(* if the matrix is large, apply Lanczos *)
		Do[
			eigs = -Eigensystem[-H, n, Method->{"Arnoldi","Criteria"->"RealPart"}];
			values = eigs[[1]];
			(* sometimes if H has degeneracies, eigenvalues are not sorted properly. We fix this by the following code *)
			If[values[[1]] > values[[2]], values = Sort[values]; eigs = SortBy[eigs\[Transpose], First]\[Transpose];];
			(* end of fixing line *)
			(* condition for stopping the loop *)
			If[T == 0,
				(* if there are no degeneracies, break the loop *)
				If[Abs[values[[-1]] - values[[-2]]] > \[Epsilon]deg, Break[];];,
			(* else if T != 0 *)
				(* if the n-th Boltzmann weight is below threshold and there are no degeneracies, break the loop *)
				If[(Abs[values[[-1]] - values[[-2]]] > \[Epsilon]deg) && (Abs[values[[-1]] - values[[1]]] > \[Epsilon]temp), Break[];];
			];
		, {n, 2, dim}];
		Return[Drop[#,-1] &/@ eigs]
	];
];
Options[Eigs] = {"Temperature" -> 0, "MinLanczosDim" -> 2, "DegeneracyThreshold" -> 10^(-9), "BoltzmannThreshold" -> 10^(-9)};

m = SparseArray[{
	Band[{1,2}]->1.,
	Band[{2,1}]->1.
},{12, 12}];
m//MatrixForm

p = SparseArray[{
	{i_,i_} :> 1.0*Mod[i,4]
},{12, 12}];


Print["Eigenstates at T=0 with several methods"];
energies = Sort@Eigenvalues[m]
Eigs[m, "MinLanczosDim" -> 20]
Eigs[m, "MinLanczosDim" -> 2]
Print["Eigenstates at T>0 with several methods"];
Eigs[m, "MinLanczosDim" -> 20, "Temperature" -> 0.1]
Eigs[m, "MinLanczosDim" -> 2, "Temperature" -> 0.1]

Print["----------------------------------------------"];

p//MatrixForm
Print["Eigenstates at T=0 with several methods"]
Sort@Eigenvalues[p]
Eigs[p, "MinLanczosDim" -> 20]
Eigs[p, "MinLanczosDim" -> 2]
Print["Eigenstates at T>0 with several methods"];
Eigs[p, "MinLanczosDim" -> 20, "Temperature" -> 0.1, "BoltzmannThreshold" -> 10^(-6)]
Eigs[p, "MinLanczosDim" -> 2, "Temperature" -> 0.1, "BoltzmannThreshold" -> 10^(-6)]


Exp[-10.*({-1.801937735804838`,-1.2469796037174674`,-0.4450418679126287`,0.4450418679126283`}-(-1.801937735804838`))]
