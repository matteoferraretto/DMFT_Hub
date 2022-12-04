(* ::Package:: *)

BeginPackage["MyLinearAlgebra`"]

Eigs::usage = "Eigs[H] returns the lowest eigenvalue(s) and the corresponding eigenvector(s) of H in the form {values, vectors}. There are several optional arguments: '
'Temperature'', ''MinLanczosDim'', ''DegeneracyThreshold'', ''BoltzmannThreshold''. 
If ''Temperature'' -> 0 (default), then the function returns only the lowest eigenstate(s) with the respective degeneracies. Two states are considered degenerate when the difference of 
their energies is below the desired ''DegeneracyThreshold'' (\!\(\*SuperscriptBox[\(10\), \(-9\)]\) is the default).
If ''Temperature'' is > 0, then all the lowest eigenstates with a Boltzmann weight above ''BoltzmannThreshold'' are returned (\!\(\*SuperscriptBox[\(10\), \(-9\)]\) is the default). 
The option ''MinLanczosDim'' establishes a threshold on the matrix dimension, above which the Lanczos method is adopted and below which full diagonalization is performed (32 by default)."

InverseElement::usage = "InverseElement[m_, {i_,j_}]"

TwoByTwoInverse::usage = "TwoByTwoInverse[A] returns the inverse of the 2x2 complex invertible matrix A. The function is listable. "


Begin["Private`"]

Print["Package MyLinearAlgebra` loaded successfully."];

(* tools to invert matrices more efficiently *)
(* get a single element of the inverse matrix *)
InverseElement[m_, {i_,j_}] := (-1)^(i+j)Det[Drop[m,{j},{i}]]/Det[m];

(* compiled and listable version of determinant specific for 2x2 matrices *)
TwoByTwoDet = Compile[{
	{A, _Complex, 2}
	},
	A[[1,1]]*A[[2,2]] - A[[1,2]]*A[[2,1]],
	CompilationTarget -> "C", RuntimeAttributes -> {Listable}, Parallelization -> True
];

(* listable inversion of 2x2 matrices *)
TwoByTwoInverse = Compile[{
	{A, _Complex, 2}
	},
	Module[
		{invA,a,b,c,d},
		a = A[[1,1]]; b = A[[1,2]]; c = A[[2,1]]; d = A[[2,2]];
		invA = (1./(a*d-b*c))*{{d,-b},{-c,a}}
	]
	, CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True
];

ThreeByThreeDet = Compile[{
	{A, _Complex, 2}
	},
	A[[1,1]]*A[[2,2]]*A[[3,3]] + A[[1,2]]*A[[2,3]]*A[[3,1]] + A[[1,3]]*A[[3,2]]*A[[2,1]] - A[[3,1]]*A[[1,3]]*A[[2,2]] - A[[3,2]]*A[[2,3]]*A[[1,1]] - A[[2,1]]*A[[1,2]]*A[[3,3]],
	CompilationTarget -> "C", RuntimeAttributes -> {Listable}, Parallelization -> True
];

(* listable inversion of 3x3 matrices *)
ThreeByThreeInverse = Compile[{
	{A, _Complex, 2}
	},
	Module[
		{trA,A2,trA2,detA,Id},
		trA = A[[1,1]]+A[[2,2]]+A[[3,3]];
		A2 = A . A;
		trA2 = A2[[1,1]]+A2[[2,2]]+A2[[3,3]];
		Id = {{1,0,0},{0,1,0},{0,0,1}};
		(* determinant according to Sarrus rule *)
		detA = A[[1,1]]*A[[2,2]]*A[[3,3]]+A[[1,2]]*A[[2,3]]*A[[3,1]]+A[[1,3]]*A[[2,1]]*A[[3,2]]
			-(A[[3,1]]*A[[2,2]]*A[[1,3]]+A[[3,2]]*A[[1,1]]*A[[2,3]]+A[[2,1]]*A[[3,3]]*A[[1,2]]);
		(* Cayley - Hamilton formula for 3x3 matrix *)
		Return[(1./detA)*(0.5*(trA^2-trA2)*Id - trA*A + A2)]
	],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True
];


(* Algorithm to compute the first element of the inverse of a tridiagonal symmetric matrix with a in the main diagonal and b in the second diagonal *)
TridiagonalInverseFirstElement = Compile[{
	{a,_Complex,1}, {b,_Complex,1}
	},
	Module[{
		n = Length[a], \[Theta] = ConstantArray[0.0, Length[a]], \[Phi] = ConstantArray[0.0, Length[a]]
	},
	(* initialize the progression \[Theta] *)
	\[Theta][[1]] = a[[1]]; 
	\[Theta][[2]] = a[[2]]*\[Theta][[1]] - b[[1]]^2;
	(* initialize the progression \[Phi] *)
	\[Phi][[n]] = a[[n]];
	\[Phi][[n-1]] = a[[n-1]]*\[Phi][[n]] - b[[n-1]]^2;
	(* compute the progression \[Theta] *)
	Do[
		\[Theta][[i]] = a[[i]]*\[Theta][[i-1]] - (b[[i-1]]^2)*\[Theta][[i-2]];
	,{i,3,n}];
	(* compute the progression \[Phi] *)
	Do[
		\[Phi][[i]] = a[[i]]*\[Phi][[i+1]] - (b[[i]]^2)*\[Phi][[i+2]];
	,{i,n-2,2,-1}];
	(* return element i=j=1 of the inverse *)
	\[Phi][[2]]/\[Theta][[n]]
	],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True
];


(* custom calculation of eigenstates *)
(* compute eigenstates *)
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
			eigs = -Eigensystem[-H, n, Method->{"Arnoldi","Criteria"->"RealPart",MaxIterations->OptionValue["MaxIterations"]}];
			values = eigs[[1]];
			(* sometimes if H has degeneracies, eigenvalues are not sorted properly. We fix this by the following code *)
			If[values[[1]] > values[[2]], values = Sort[values]; eigs = SortBy[eigs\[Transpose], First]\[Transpose];];
			(* end of fixing line *)
			(* condition for stopping the loop *)
			If[T == 0,
				(* if there are no degeneracies, break the loop *)
				If[Abs[values[[-1]] - values[[1]]] > \[Epsilon]deg, Break[];];,
			(* else if T != 0 *)
				(* if the n-th Boltzmann weight is below threshold and there are no degeneracies, break the loop *)
				If[(Abs[values[[-1]] - values[[-2]]] > \[Epsilon]deg) && (Abs[values[[-1]] - values[[1]]] > \[Epsilon]temp), Break[];];
			];
		, {n, OptionValue["MinEigenvalues"], dim}];
		Return[Drop[#,-1] &/@ eigs]
	];
];
Options[Eigs] = {"Temperature" -> 0, "MinLanczosDim" -> 32, "DegeneracyThreshold" -> 10^(-9), "BoltzmannThreshold" -> 10^(-9), "MaxIterations" -> 1000, "MinEigenvalues" -> 10};


(*                       LANCZOS                       *)
Lanczos[H_, StartingVector_, OptionsPattern[]] := Module[{
	miniter = OptionValue[MinIter],
	maxiter = Min[OptionValue[MaxIter], Length[H]-1],
	\[Epsilon] = OptionValue[ConvergenceThreshold],
	shift = OptionValue[Shift],
	a,b,a0,b1,v,w,HKrilov,E0old,E0new,nfinal
	},
	(* initialize array of a_n :  a[1]=a_0 , a[1+n]=a_n *)
	a = ConstantArray[0, maxiter+1];
	(*initialize array of b_n : b[n]=b_n *)
	b = ConstantArray[0, maxiter];
	(* initialize starting vector *)
	v = Normalize[StartingVector];
	(* initialize a new vector w = Hv *)
	w = H . v;
	a0 = (Conjugate[v]) . w;
	w = w - a0*v;
	b1 = Norm[w];
	a[[1]] = a0; b[[1]] = b1;
	E0old = a0;
	nfinal = maxiter;
	Do[
		If[b[[n]] < \[Epsilon] && n >= miniter,
			nfinal = n; Break[];
		];
		w = w/b[[n]];(* w=Subscript[v, n] *)
		v = -b[[n]]*v;(* v=-Subscript[b, n]Subscript[v, n-1]*)
		{v,w} = {w,v};
		w = w + H . v;(* w=Subscript[Hv, n]-Subscript[b, n]Subscript[v, n-1] *)
		a[[n+1]] = (Conjugate@v) . w; (*Subscript[a, n] = Subscript[v, n]Subscript[Hv, n]-Subscript[b, n]Subscript[v, n]Subscript[v, n-1]*)
		w = w - a[[n+1]] * v;
		If[n < maxiter,
			b[[n+1]] = Norm[w];
		];
		(*Build hamiltonian in the Krylov subspace*)
		HKrilov = SparseArray[{
			Band[{1,2}] -> Take[b, n],
			Band[{2,1}] -> Take[b, n],
			Band[{1,1}] -> Take[a, n+1]
		},{n+1, n+1}];
		E0new = If[shift == 0,
			Min @ Eigenvalues[HKrilov, Method -> "Banded"],
		(*else*)	
			Eigenvalues[HKrilov-DiagonalMatrix[ConstantArray[shift,1+n]]][[1]]+shift
		];
		If[Abs[E0new - E0old] < \[Epsilon] && n > miniter,
			nfinal = n; Break[];
		];
		E0old = E0new;
	, {n, 1, maxiter}];
	a = Take[a, nfinal+1];
	b = Take[b, nfinal];
	{E0new,a,b}
];
Options[Lanczos] = {ConvergenceThreshold -> 1.0*10^(-8), MinIter -> 1, MaxIter -> 2000, Shift -> 0};

(* avoid eigenvalue computation inside, and performed a fixed number of iterations. 
WARNING: it is faster than Lanczos[] for a fixed number of iterations, yet it can be overall slower since it never stops earlier *)
LanczosCore[H_, StartingVector_, OptionsPattern[]] := Module[{
	maxiter = Min[OptionValue[MaxIter], Length[H]-1],
	\[Epsilon] = OptionValue[ConvergenceThreshold],
	a,b,a0,b1,v,w,nfinal,HKrylov,E0
	},
	(* initialize array of a_n :  a[1]=a_0 , a[1+n]=a_n *)
	a = ConstantArray[0.0, maxiter+1];
	(*initialize array of b_n : b[n]=b_n *)
	b = ConstantArray[0.0, maxiter];
	(* initialize starting vector *)
	v = Normalize[StartingVector];
	(* initialize a new vector w = Hv *)
	w = H . v;
	a0 = (Conjugate[v]) . w;
	w = w - a0*v;
	b1 = Norm[w];
	a[[1]] = a0; b[[1]] = b1;
	nfinal = maxiter;
	(* start the loop *)
	Do[
		If[b[[n]] < \[Epsilon],
			nfinal = n;
			Print["Lanczos terminated after "<>ToString[n]<>" out of "<>ToString[maxiter]<>" iterations."];
			Break[];
		];
		w = w/b[[n]];(* w=Subscript[v, n] *)
		v = -b[[n]]*v;(* v=-Subscript[b, n]Subscript[v, n-1]*)
		{v,w} = {w,v};
		w = w + H . v;(* w=Subscript[Hv, n]-Subscript[b, n]Subscript[v, n-1] *)
		a[[n+1]] = (Conjugate@v) . w; (*Subscript[a, n] = Subscript[v, n]Subscript[Hv, n]-Subscript[b, n]Subscript[v, n]Subscript[v, n-1]*)
		w = w - a[[n+1]] * v;
		If[n < maxiter,
			b[[n+1]] = Norm[w];
		];
	, {n, 1, maxiter}];
	a = a[[;;nfinal+1]];
	b = b[[;;nfinal]];
	HKrylov = SparseArray[{
			Band[{1,2}] -> b,
			Band[{2,1}] -> b,
			Band[{1,1}] -> a
		}, {nfinal+1, nfinal+1}];
	E0 = Min[Eigenvalues[HKrylov, Method -> "Banded"]];
	{E0, a, b}
];
Options[LanczosCore] = {ConvergenceThreshold -> 1.0*10^(-8), MaxIter -> 2000};

End[]

EndPackage[]
