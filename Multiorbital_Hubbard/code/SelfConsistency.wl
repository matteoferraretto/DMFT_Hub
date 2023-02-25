(* ::Package:: *)

BeginPackage["SelfConsistency`", {"MyLinearAlgebra`"}]


(* Self consistency *)
WeissField::usage = "WeissField[L_, f_, Norb_, \[Mu]_, symbols_, z_, EdMode]"

WeissFieldNumeric::usage = "WeissFieldNumeric[DBethe, \[Mu], LocalG, LocalGold, \[CapitalSigma], \[CapitalSigma]old, zlist, EdMode, OptionsPattern] gives a numeric evaluation of the Weiss field at the n-th iteration.
The calculation is based on the knowledge of the n-th and (n-1)-th local Green function and the n-th and (n-1)-th self energy, that are mixed up when setting the option Mix -> ...
The calculation depends on the lattice type: it takes advantage of the simplified algebra when Lattice -> ''Bethe'' and LatticeDimension -> Infinity; while it is Gloc^-1 + \[CapitalSigma] in general.
The list of Matsubara frequencies zlist must be provided, along with the effective chemical potential. The parameter DBethe represents the half-bandwidth in case of Bethe lattice, while
it is ignored when the lattice type is not ''Bethe''. "

SelfConsistency::usage = "SelfConsistency[Weiss, symbols, z, IndependentParameters, WeissNumeric, zlist, EdMode]" 

DMFTError::usage = "DMFTError[Xnew, Xold, EdMode]"


Begin["Private`"]

Print["Package SelfConsistency` loaded successfully."];

(* ANALYTIC EVALUATION OF NONINTERACTING IMPURITY GREEN FUNCTION *)
(* symbolic non-interacting Green function obtained by inverting the impurity-bath Hamiltonian *)
WeissField[L_, f_, Norb_, \[Mu]_, symbols_, z_, EdMode_] := Module[
	{e, V, \[CapitalDelta], \[CapitalXi], H0, H, Vmat},
	Which[
		EdMode == "Normal", 
		e = symbols[[1;;L-1]];
		V = symbols[[L;;2(L-1)]];
		z + \[Mu] - Sum[V[[k]]^2/(z - e[[k]]), {k, L-1}],
	(* ---------------------------------------- *)
		EdMode == "Superc", 
		e = symbols[[1;;L-1]];
		V = symbols[[L;;2(L-1)]];
		\[CapitalDelta] = symbols[[2L-1;;3(L-1)]];
		(* initialize stuff *)
		H = ConstantArray[0, {L-1, 2, 2}];
		Vmat = H;
		(* define the hamiltonian blocks *)
		H0 = DiagonalMatrix[{-\[Mu], + \[Mu]}];
		(H[[#]] = e[[#]]*PauliMatrix[3] + \[CapitalDelta][[#]]*PauliMatrix[1]) &/@ Range[L-1];
		(Vmat[[#]] = V[[#]]*PauliMatrix[3]) &/@ Range[L-1];
		(* compute Weiss field *)
		z*IdentityMatrix[2] - H0 - Sum[Vmat[[k]] . Inverse[z*IdentityMatrix[2] - H[[k]]] . Vmat[[k]] , {k, L-1}],
	(* ---------------------------------------- *)
		EdMode == "Raman",
		e = Partition[ symbols[[ ;; (L-1)*f*(f+1)/2]], f*(f+1)/2];
		V = Partition[ symbols[[-(L-1)*f*(f+1)/2 ;; ]], f*(f+1)/2];
		(* initialize stuff *)
		H = ConstantArray[0, {L-1, f, f}];
		Vmat = H;
		(* define the hamiltonian blocks *)
		H0 = \[Mu] * IdentityMatrix[f];
		(H[[#]] = SymmetricMatrixFromArray[e[[#]], f]) &/@ Range[L-1];
		(Vmat[[#]] = SymmetricMatrixFromArray[V[[#]], f]) &/@ Range[L-1];
		(* compute Weiss field *)
		z*IdentityMatrix[f] - H0 - Sum[Vmat[[k]] . Inverse[z*IdentityMatrix[f] - H[[k]]] . Vmat[[k]] , {k, L-1}],
	(* ---------------------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		(* split symbols according to their meaning, e, V ... *)
		e = symbols[[1 ;; Norb*(L-1)]];
		V = symbols[[Norb*(L-1)+1 ;; 2*Norb*(L-1)]];
		\[CapitalDelta] = symbols[[2*Norb*(L-1)+1 ;; 3*Norb*(L-1)]];
		\[CapitalXi] = symbols[[3*Norb*(L-1)+1 ;; (3*Norb+1)*(L-1)]];
		(* separate into the various orbital parts *)
		e = Partition[e, L-1]; V = Partition[V, L-1]; \[CapitalDelta] = Partition[\[CapitalDelta], L-1];
		(* get Hamiltonian blocks *)
		(* the chemical potential is now a list of Norb numbers given by \[Mu] + \[Delta][[orb]] *)
		(* Only compatible with Mathematica 13.1 *)
		(* H0 = BlockDiagonalMatrix[
			Table[-\[Mu][[orb]] * PauliMatrix[3], {orb, Norb}]
		]; *)
		(* compatible with Mathematica 13.0 *)
		H0 = DiagonalMatrix[
			Flatten[
				Table[-\[Mu][[orb]]*{1, -1}, {orb, Norb}]
			]
		];
		H = SparseArray[{}, {L-1, Norb, Norb}];
		Vmat = SparseArray[{}, {L-1, Norb, Norb}];
		Do[
			If[orb1 == orb2,
				H[[#]][[orb1, orb1]] = e[[orb1]][[#]] * PauliMatrix[3] + \[CapitalDelta][[orb1]][[#]] * PauliMatrix[1];
				Vmat[[#]][[orb1, orb1]] = V[[orb1]][[#]] * PauliMatrix[3];
			];
		(* else, if orb1 != orb2 *)
			If[orb2 > orb1,
				H[[#]][[orb1, orb2]] = \[CapitalXi][[#]] * PauliMatrix[1]
			];
			If[orb2 < orb1,
				H[[#]][[orb1, orb2]] = \[CapitalXi][[#]] * PauliMatrix[1]
			];
		, {orb1, Norb}, {orb2, Norb}] &/@ Range[L-1];
		(* give a matrix form to the block form above *)
		Do[
			H[[k]] = ArrayFlatten[H[[k]]];
			Vmat[[k]] = ArrayFlatten[Vmat[[k]]];
		, {k, L-1}];
		(* compute Weiss field *)
		z * IdentityMatrix[2Norb] - H0 - Sum[Vmat[[k]] . (Inverse[z * IdentityMatrix[2Norb] - H[[k]]] . Vmat[[k]]) , {k, L-1}]
	]
];

(* numeric evaluation of the Weiss field, including the mixing with previous iteration *)
WeissFieldNumeric[DBethe_, \[Mu]_, LocalG_, LocalGold_, \[CapitalSigma]_, \[CapitalSigma]old_, zlist_, EdMode_, OptionsPattern[]] := Module[
	{\[Alpha] = OptionValue[Mix], Lattice = OptionValue[Lattice], d = OptionValue[LatticeDimension], \[Sigma]3 = PauliMatrix[3], LocalGeff, Weff},
	Which[
		EdMode == "Normal" && Lattice == "Bethe" && d == Infinity,
		(* Mix up the local Green function *)
		LocalGeff = \[Alpha] * LocalGold + (1.0 - \[Alpha]) * LocalG;
		Weff = zlist + \[Mu] - (DBethe^2/4.)*LocalGeff;,
	(* ------------------------------------ *)
		EdMode == "Superc" && Lattice == "Bethe" && d == Infinity,
		(* Mix up the local Green function *)
		LocalGeff = \[Alpha] * LocalGold + (1.0 - \[Alpha]) * LocalG;
		Weff = (#*IdentityMatrix[2] + \[Mu]*\[Sigma]3) &/@ zlist - (DBethe^2/4.) * Map[Dot[\[Sigma]3, #, \[Sigma]3]&, LocalGeff];,
	(* ------------------------------------ *)
		EdMode == "Normal" && Lattice != "Bethe",
		(* Mix up the effective Weiss field *)
		If[\[Alpha] == 0.0,
			Weff = 1./LocalG + \[CapitalSigma];,
		(* else, if mixing is active *)
			Weff = \[Alpha] * (1./LocalGold + \[CapitalSigma]old) + (1.0 - \[Alpha]) * (1./LocalG + \[CapitalSigma]);
		];,
	(* ------------------------------------ *)
		EdMode == "Superc" && Lattice != "Bethe",
		(* Mix up the effective Weiss field *)
		If[\[Alpha] == 0.0,
			Weff = TwoByTwoInverse[LocalG] + \[CapitalSigma];,
		(* else, if mixing is active *)
			Weff = \[Alpha] * (TwoByTwoInverse[LocalGold] + \[CapitalSigma]old) + (1.0 - \[Alpha]) * (TwoByTwoInverse[LocalG] + \[CapitalSigma]);
		];,
	(* ------------------------------------ *)
		EdMode == "Raman",
		If[\[Alpha] == 0.0,
			Weff = Inverse[LocalG] + \[CapitalSigma];,
		(* else, if mixing is active *)
			Weff = \[Alpha] * (Inverse[LocalGold] + \[CapitalSigma]old) + (1.0 - \[Alpha]) * (Inverse[LocalG] + \[CapitalSigma]);
		],
	(* ------------------------------------ *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		(* Mix up the effective Weiss field *)
		If[\[Alpha] == 0.0, 
			Weff = (Inverse /@ LocalG) + \[CapitalSigma];, (* if no mixing, avoid inverting the old local G *)
		(* else, if mixing is active *)
			Weff = \[Alpha] * ((Inverse /@ LocalGold) + \[CapitalSigma]old) + (1.0 - \[Alpha]) * ((Inverse /@ LocalG) + \[CapitalSigma]);
		];
	];
	Weff
];
Options[WeissFieldNumeric] = {Mix -> 0.0, Lattice -> "Bethe", LatticeDimension -> Infinity};

(* new self cons *)
SelfConsistency[Weiss_, symbols_, z_, IndependentParameters_, WeissNumeric_, zlist_, EdMode_, OptionsPattern[{SelfConsistencyNew, FindMinimum}]] := Module[
	{weight = OptionValue[FitWeight], \[Chi], residue, newparameters},
	Which[
		EdMode == "Normal",
		(* distance function *)
		\[Chi][symbols] = Mean[Abs[ weight * (
			(Weiss /.{z -> #} &/@ zlist) - WeissNumeric
		)]^2],
	(* ------------------------------------ *)
		EdMode == "Superc",
		(* distance function *)
		\[Chi][symbols] = Mean @ First @
			Mean[Abs[ weight * (
				(Weiss/.{z -> #} &/@ zlist) - WeissNumeric
			)]^2],
	(* ------------------------------------ *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		(* distance function *)
		\[Chi][symbols] = Mean[
			Total[#, 2] &/@ (
				Abs[
					weight * (
					(UpperTriangularize[Weiss/.{z -> #}] &/@ zlist)
					- (UpperTriangularize[#] &/@ WeissNumeric)
				)]^2 
			)];
	];
	(* perform the fit to find the new bath parameters *)
	Which[
		OptionValue[Minimum] == "Local",
		{residue, newparameters} =
			FindMinimum[
				\[Chi][symbols],
				{symbols, IndependentParameters}\[Transpose],
				Method -> OptionValue[Method],
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
Options[SelfConsistency] = { 
	NumberOfFrequencies -> 2000, 
	Minimum -> "Local", 
	Method -> "ConjugateGradient", 
	FitWeight -> ConstantArray[1., 2000]
};


(* Perform minimization according to the self consistency condition *)
SelfConsistencyOld[DBethe_, \[Mu]_, Weiss_, symbols_, z_, IndependentParameters_, LocalG_, LocalGold_, \[CapitalSigma]_, \[CapitalSigma]old_, zlist_, EdMode_, OptionsPattern[{SelfConsistency, FindMinimum}]] := 
Module[
	{Lattice = OptionValue[Lattice],
	weight = OptionValue[FitWeight],
	\[Alpha] = OptionValue[Mix],
	LocalGeff, Weff,
	residue, newparameters, \[Chi], \[Sigma]3 = PauliMatrix[3]},
	(* define the target function to minimize depending on the Lattice and EdMode (if Lattice = Bethe there is a shortcut) *)
	Which[
		EdMode == "Normal" && Lattice == "Bethe",
		(* Mix up the local Green function *)
		LocalGeff = \[Alpha] * LocalGold + (1.0 - \[Alpha]) * LocalG;
		(* distance function *)
		\[Chi][symbols] = Mean[Abs[ weight * (
			((Weiss - z - \[Mu])/.{z -> #} &/@ zlist) + (DBethe^2/4.)*LocalGeff
		)]^2],
	(* ------------------------------------ *)
		EdMode == "Superc" && Lattice == "Bethe",
		(* Mix up the local Green function *)
		LocalGeff = \[Alpha] * LocalGold + (1.0 - \[Alpha]) * LocalG;
		(* distance function *)
		\[Chi][symbols] = Mean @ First @
			Mean[Abs[ weight *
				((Weiss - z - \[Mu]*\[Sigma]3)/.{z -> #} &/@ zlist) 
				+ (DBethe^2/4.) * Map[Dot[\[Sigma]3, #, \[Sigma]3]&, LocalGeff]
			]^2],
	(* ------------------------------------ *)
		EdMode == "Normal" && Lattice != "Bethe",
		(* Mix up the effective Weiss field *)
		If[\[Alpha] == 0.0,
			Weff = 1./LocalG + \[CapitalSigma],
		(* else, if mixing is active *)
			Weff = \[Alpha] * (1./LocalGold + \[CapitalSigma]old) + (1.0 - \[Alpha]) * (1./LocalG + \[CapitalSigma])
		];
		(* distance function *)
		\[Chi][symbols] = Mean[Abs[ weight * (
			(Weiss/.{z -> #} &/@ zlist) - Weff
		)]^2],
	(* ------------------------------------ *)
		EdMode == "Superc" && Lattice != "Bethe",
		(* Mix up the effective Weiss field *)
		If[\[Alpha] == 0.0,
			Weff = TwoByTwoInverse[LocalG] + \[CapitalSigma],
		(* else, if mixing is active *)
			Weff = \[Alpha] * (TwoByTwoInverse[LocalGold] + \[CapitalSigma]old) + (1.0 - \[Alpha]) * (TwoByTwoInverse[LocalG] + \[CapitalSigma])
		];
		(* distance function *)
		\[Chi][symbols] = Mean @ First @
			Mean[Abs[
				(Weiss/.{z -> #}& /@ zlist) - Weff
			]^2],
	(* ------------------------------------ *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		(* Mix up the effective Weiss field *)
		If[\[Alpha] == 0.0, 
			Weff = (Inverse /@ LocalG) + \[CapitalSigma], (* if no mixing, avoid inverting the old local G *)
		(* else, if mixing is active *)
			Weff = \[Alpha] * ((Inverse /@ LocalGold) + \[CapitalSigma]old) + (1.0 - \[Alpha]) * ((Inverse /@ LocalG) + \[CapitalSigma])
		];
		(* distance function *)
		\[Chi][symbols] = Mean[
			Total[#, 2] &/@ (
				Abs[
					weight * (
					(UpperTriangularize[Weiss/.{z -> #}] &/@ zlist) 
					- (UpperTriangularize[#] &/@ Weff)
				)]^2 
			)];
	];
	(* perform the fit to find the new bath parameters *)
	Which[
		OptionValue[Minimum] == "Local",
		{residue, newparameters} =
			FindMinimum[
				\[Chi][symbols],
				{symbols, IndependentParameters}\[Transpose],
				Method -> OptionValue[Method],
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
Options[SelfConsistencyOld] = {
	Mix -> 0.0,
	Lattice -> "Bethe", 
	LatticeDimension -> 2, 
	NumberOfFrequencies -> 2000, 
	Minimum -> "Local", 
	Method -> "ConjugateGradient", 
	FitWeight -> ConstantArray[1., 2000]
};


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
(* --------------------------------- *)
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
		}],
(* ----------------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		error = {};
		Do[
			If[EvenQ[\[Alpha]] && \[Alpha] == \[Beta], Continue[]; ]; (* skip elements 22, 44, ... because they are related to 11, 33, ... *)
			If[
				Max[Total[Abs[Xnew[[All, \[Alpha], \[Beta]]]]], Total[Abs[Xold[[All, \[Alpha], \[Beta]]]]]] < 0.1,
				Continue[]; 
			]; (* skip elements when they are zero *)
			Print["\[Alpha]=",\[Alpha]," \[Beta]=",\[Beta]];
			AppendTo[error, Total[
					Abs[Xnew[[All, \[Alpha], \[Beta]]] - Xold[[All, \[Alpha], \[Beta]]]]
				]/Max[
					Total[Abs[Xnew[[All, \[Alpha], \[Beta]]]]], Total[Abs[Xold[[All, \[Alpha], \[Beta]]]]]		
				]
			];
		, {\[Alpha], Length[Xnew[[1]]]}, {\[Beta], \[Alpha], Length[Xnew[[1]]]}];
		error = Mean[error];
		(*
		error = Mean[
			Flatten @ Table[
				(* else *)
					Total[
						Abs[Xnew[[All, \[Alpha], \[Beta]]] - Xold[[All, \[Alpha], \[Beta]]]]
					]/Max[
						Total[Abs[Xnew[[All, \[Alpha], \[Beta]]]]], Total[Abs[Xold[[All, \[Alpha], \[Beta]]]]]		
					]	
			, {\[Alpha], Length[Xnew[[1]]]}, {\[Beta], \[Alpha], Length[Xnew[[1]]]}]
		]*)
	];
	error
];

End[]

EndPackage[]
