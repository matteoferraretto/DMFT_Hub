(* ::Package:: *)

Eigs[H_, OptionsPattern[]] := Module[
	{values, vectors, eigs, dim = Length[H], T = OptionValue["Temperature"], \[Epsilon]deg = OptionValue["DegeneracyThreshold"]},
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
						values, x_/;Abs[x - Min[values]] < \[Epsilon]deg
					]
				]&/@eigs	
			],
		(* else *)
			Return["todo"]
		],
	(* else *)
		(* if the matrix is large, apply Lanczos *)
		Do[
			eigs = -Eigensystem[-H, n, Method->{"Arnoldi","Criteria"->"RealPart"}];
			values = eigs[[1]];
			(* sometimes if H has degeneracies, eigenvalues are not sorted properly. We fix this by the following code *)
				If[values[[1]] > values[[2]], values = Sort[values]; eigs = SortBy[eigs\[Transpose], First]\[Transpose];];
			(* end of fixing line *)
			(* if there are no degeneracies, break the loop *)
			If[Abs[values[[-1]] - values[[-2]]] > \[Epsilon]deg, Break[];];
		, {n, 2, dim}];
		Return[Drop[#,-1] &/@ eigs];
	];
];
Options[Eigs] = {"Temperature" -> 0, "MinLanczosDim" -> 2, "DegeneracyThreshold" -> 10^(-9)};

m = SparseArray[{
	Band[{1,2}]->1.,
	Band[{2,1}]->1.
},{6,6}];
m//MatrixForm

p = SparseArray[{
	{i_,i_} :> 1.0*Mod[i,3]
},{8,8}];
p//MatrixForm

Sort@Eigenvalues[m]
Eigs[m, "MinLanczosDim" -> 20]
Eigs[m, "MinLanczosDim" -> 2]

Print["----------------------------------------------"];

Sort@Eigenvalues[p]
Eigs[p, "MinLanczosDim" -> 20]
Eigs[p, "MinLanczosDim" -> 2]



