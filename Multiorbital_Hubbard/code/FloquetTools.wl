(* ::Package:: *)

BeginPackage["FloquetTools`"]



Begin["Private`"]

Print["Package FloquetTools` loaded successfully."];

(* time contour in the complex plane *)
Contour = Compile[{
	{tmin,_Real}, {tmax,_Real}, {\[Beta],_Real}, {Nt,_Integer}
	},
	With[{dt = 1.0*(tmax - tmin)/Nt},
		Join[
			Table[tmin + (j - 1)*dt + 0.0*I, {j, 1, Nt}],
			Table[tmax - (j - 1 - Nt)*dt + 0.0*I, {j, Nt+1, 2Nt}],
			Table[tmin - (\[Beta]/Nt)*I*(j - 1 - 2Nt), {j, 2Nt+1, 3Nt}]
		]
	]
];

(* plot the contour *)
ShowContour[contour_] := With[
	{Nt = Length[contour]/3, tmin = contour[[1]], \[Eta] = 0.05},
	Print @ ListPlot[{
		{contour[[;;Nt]], ConstantArray[\[Eta], Nt]}\[Transpose],
		{contour[[Nt+1;;2Nt]], ConstantArray[-\[Eta], Nt]}\[Transpose],
		{ConstantArray[Re[tmin], Nt], Im[contour[[2Nt+1;;3Nt]]]}\[Transpose]
	}]
];

(* weights along the contour for computing integrals *)
ContourWeights = Compile[{
		{contour, _Complex, 1}
	},
	With[
		{Nt = Length[contour]/3, dt = contour[[2]]-contour[[1]], d\[Tau] = contour[[-2]]-contour[[-1]]},
		Join[
			Table[dt, {n, Nt}],
			Table[-dt, {n, Nt}],
			Table[-d\[Tau], {n, Nt}]
		]
	]
];

(* integrate a function over the contour *)
ContourIntegrate = Compile[{
		{contour,_Complex,1}, {fvalues,_Complex,1}
	},
	With[
		{W = ContourWeights[contour]},
		W . fvalues
	]
];

(* define the contour-delta function \[Delta](t, t') as an off-diagonal operator, which is diagonal when dt -> 0 *)
(* this is the definition when \[Delta] is to the left in the integrand (integral over j) *)
ContourDeltaLeft = Compile[{
		{contour, _Complex, 1}, {i, _Integer}, {j, _Integer}
	},
	With[{W = ContourWeights[contour]},
		If[j == Length[contour] && i == 1, Return[-1./W[[i]]] ];
		If[
			i == j + 1, 1./W[[i]],
		(* else *)	
			0
		]
	], CompilationTarget -> "C", RuntimeAttributes -> {Listable}, Parallelization -> True
];

(* this is the definition when \[Delta] is to the left in the integrand *)
ContourDeltaRight = Compile[{
		{contour, _Complex, 1}, {i, _Integer}, {j, _Integer}
	},
	With[{W = ContourWeights[contour]},
		If[i == 1 && j == Length[contour], Return[-1./W[[j]]] ];
		If[
			i - 1 == j, 1./W[[j]],
		(* else *)	
			0
		]
	], CompilationTarget -> "C", RuntimeAttributes -> {Listable}, Parallelization -> True
];

(* returns a matrix representation of the operator (i d/dt + \[Mu])\[Delta](t_j, t_k) *)
ContourDeltaDerivative[contour_, \[Mu]_] := With[
	{Nt = Length[contour]/3, dt = contour[[2]]-contour[[1]], d\[Tau] = contour[[-2]]-contour[[-1]], W = ContourWeights[contour]},
	I * SparseArray[{
		Band[{1,1}] -> 1./W^2,
		Band[{2,1}] -> Join[
			ConstantArray[-1. - I*dt*\[Mu], Nt-1],
			ConstantArray[-1. + I*dt*\[Mu], Nt],
			ConstantArray[-1. - Im[d\[Tau]]*\[Mu], Nt]
		] * (1./W[[2;;]]) * (1./W[[;;-2]]),
		{1, 3Nt} -> (1. + I*dt*\[Mu])/(W[[1]]*W[[-1]])
	}
	, {3Nt, 3Nt}]
];

(* M_jk matrix: it is defined via the relation (id_t + \[Mu])\[Delta](t_j, t_k) = i M_jk/(W_j W_k) *)
ContourM[contour_, \[Mu]_] := With[
	{Nt = Length[contour]/3, dt = contour[[2]]-contour[[1]], d\[Tau] = contour[[-2]]-contour[[-1]], W = ContourWeights[contour]},
	SparseArray[{
		Band[{1,1}] -> 1.,
		Band[{2,1}] -> Join[
			ConstantArray[-1. - I*dt*\[Mu], Nt-1],
			ConstantArray[-1. + I*dt*\[Mu], Nt],
			ConstantArray[-1. - Im[d\[Tau]]*\[Mu], Nt]
		],
		{1, 3Nt} -> (1. + I*dt*\[Mu])
	}
	, {3Nt, 3Nt}]
];





End[]

EndPackage[]
