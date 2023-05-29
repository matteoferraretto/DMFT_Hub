(* ::Package:: *)

BeginPackage["Lattices`", {"MyLinearAlgebra`"}]


DoSBethe::usage = "DoSBethe[\[Epsilon], DBethe] returns the density of states for the infinite dimensional Bethe lattice with half bandwidth DBethe and energy \[Epsilon]"

BrillouinZone::usage = "BrillouinZone[LE, LatticeDim, LatticeType] returns a list of LE^LatticeDim points in the first Brillouin zone for a given lattice, 
for example LatticeType=''Hypercubic''. LatticeDim represents the dimensionality of such lattice, hence the dimension of each point in the Brillouin zone. "

MagneticBrillouinZone::usage = "MagneticBrillouinZone[LE, LatticeDim, LatticeType] returns a list of (LE/2)^LatticeDim points in the first Magnetic Brillouin zone 
(i.e. reduced Brillouin zone when the lattice Hamiltonian is written in the basis of two sublattices of the original lattice).  "

GetLatticeEnergies::usage = "GetLatticeEnergies[HalfBandwidths, \[Delta], LatticeType, LatticeDim, NumberOfPoints] returns a list with two objects:
{LatticeEnergies, LatticeWeights}. LatticeWeights is a list of weights used when performing sums over the momenta/energies; if summing over momenta this is simply a uniform list that
is normalized to 1; if summing over energies this is a sampling of the density of states. LatticeEnergies is a nested list with dimension {NumberOfPoints, Norb, Norb}. 
The first index labels the momentum value; for every momentum k we have a matrix \[Epsilon]_orb1,orb2(k). Here \[Delta] is the crystal field splitting."

GetLatticeEnergiesRaman::usage = "GetLatticeEnergiesRaman[HalfBandwidths, \[Delta], M, \[Gamma], u, LatticeType, LatticeDim, NumberOfPoints] returns a list with two objects:
{LatticeEnergies, LatticeWeights}. LatticeWeights is a list of weights used when performing sums over the momenta/energies; if summing over momenta this is simply a uniform list that
is normalized to 1; if summing over energies this is a sampling of the density of states. LatticeEnergies is a nested list with dimension {NumberOfPoints, Norb, Norb, f, f}. 
The first index labels the momentum value; for every momentum k we have a matrix \[Epsilon]_orb1,orb2,\[Sigma],\[Rho](k) written in the basis where the Raman matrix M is diagonal. Here \[Gamma] and u represent the
magnitude and direction of the gauge field, while \[Delta] is the crystal field splitting. "

GetLatticeEnergiesRamanSublattices::usage = "GetLatticeEnergiesRamanSublattices[HalfBandwidths, \[Delta], M, \[Gamma], LatticeType, LatticeDim, NumberOfPoints]"

DispersionHypercubic::usage = "."

DispersionHypercubicRaman::usage = "DispersionHypercubicRaman[k, t, f, \[Gamma], u] returns a fxf matrix that represents the energy associated to the given point k in the 1st Brillouin zone
in presence of a gauge field of magnitude \[Gamma] along the direction specified by the unit vector u. The gauge field shifts anti-symmetrically the dispersion relations for different spin states. The input
t is the nearest neighbor hopping energy scale. Notice that there is no Raman tunneling in these matrices."

HighSymmetryPath::usage = "."

LocalGreenFunction::usage = "LocalGreenFunction[LatticeEnergies, weights, \[Mu], \[CapitalSigma], zlist, EdMode]"


Begin["Private`"]

Print["Package Lattices` loaded successfully."];

(* density of states of the infinite dimensional Bethe lattice *)
DoSBethe = Compile[
	{{\[Epsilon], _Real}, {DBethe, _Real}},
	(2./(Pi*DBethe^2))*Sqrt[DBethe^2 - \[Epsilon]^2], 
	CompilationTarget->"C", RuntimeAttributes->{Listable}
];

(* dispersion relation for a d-dimensional hypercubic lattice *)
DispersionHypercubic = Compile[
	{{k,_Real,1}, {t,_Real}},
	-2.*t*Sum[Cos[ka], {ka, k}],
	CompilationTarget -> "C", RuntimeAttributes -> {Listable}
];

(* dispersion relation for a d-dimensional hypercubic lattice in presence of Raman field with gauge field *)
DispersionHypercubicRaman = Compile[
	{{k,_Real,1}, {t,_Real}, {M,_Real,2}, {\[Gamma],_Real,1}},
	-2.*t*DiagonalMatrix[
		Table[
			Sum[Cos[k[[a]] - \[Sigma] * \[Gamma][[a]]], {a, Length[k]}]
		, {\[Sigma], -(Length[M]-1)/2, (Length[M]-1)/2}]
	] + M,
	CompilationTarget -> "C", RuntimeAttributes -> {Listable}
];


(* returns a list of k points (d-dimensional vectors) in the 1st BZ *)
BrillouinZone[LE_, LatticeDim_, LatticeType_] := With[
	{dk = 2.Pi/LE, d = LatticeDim},
	Which[
		LatticeType == "Hypercubic",
		Tuples[ Table[k, {k, -1.*Pi+dk, 1.*Pi, dk}], d],
	(* ----------------------------------------------- *)
		LatticeType != "Hypercubic",
		Print["Brillouin zone for the specified lattice is not supported."];
	]
];

(* Magnetic Brillouin zone *)
MagneticBrillouinZone[LE_, LatticeDim_, LatticeType_] := With[
	{BZ = BrillouinZone[LE, LatticeDim, LatticeType], d = LatticeDim},
	If[
		LatticeType == "Hypercubic",
		Which[
			d == 1,
			Select[BZ, (#[[1]] <= Pi/2 && #[[1]] > -Pi/2)&],
		(* ----------------------------------- *)
			d == 2,
			Select[BZ, (#[[2]]<=Pi-Abs[#[[1]]] && #[[2]]>-Pi+Abs[#[[1]]])&]
		],
		(* else *)
		Print["Magnetic Brillouin zone for the specified lattice is not supported."]
	]
];


(* "high symmetry path" for the Bethe lattice. This is an abuse of notation: this just returns the indexes of all energies *)
HighSymmetryPathBethe[LatticePoints_] := Range[LatticePoints];

(* high symmetry path \[CapitalGamma]XM\[CapitalGamma] for the square lattice *)
HighSymmetryPathSquare[LatticePoints_] := Module[
	{LE = Sqrt[LatticePoints], \[CapitalGamma], X, M},
	If[OddQ[LE], 
		Print["Error: an odd number of points per lattice dimension does not describe the high symmetry path. "]; 
	];
	\[CapitalGamma] = LE*(LE/2-1) + LE/2;
	X = \[CapitalGamma] + (LE^2)/2;
	M = LatticePoints;
	Join[
		Table[\[CapitalGamma] + n*LE, {n, 0, LE/2}],
		Table[X + n, {n, 1, LE/2}],
		Table[M - n*(LE+1), {n, 1, LE/2}]
	]
];

HighSymmetryPath[LatticePoints_, LatticeType_, LatticeDim_] := Which[
	LatticeType == "Bethe", HighSymmetryPathBethe[LatticePoints],
	LatticeType == "Hypercubic" && LatticeDim == 1, HighSymmetryPathBethe[LatticePoints], (* if d=1, the high symmetry path is just like in Bethe lattice *)
	LatticeType == "Hypercubic" && LatticeDim == 2, HighSymmetryPathSquare[LatticePoints],
	True, Print["Not supported."];
];

(* split the square lattice into two sublattices *)
(*SplitSquareLatticeIntoSublattices[BZ_] := Module[
	{NumberOfPoints, LE, d, row, col, indexes, BZA, BZB},
	{NumberOfPoints, d} = Dimensions[BZ];
	LE = NumberOfPoints^(1/d);
	indexes = Range[NumberOfPoints];
	row = Table[Mod[indexes[[i]]-1, LE], {i, NumberOfPoints}];
	col = Table[Quotient[indexes[[i]]-1, LE], {i, NumberOfPoints}];
	BZA = Pick[BZ, EvenQ[row + col]];
	BZB = Pick[BZ, OddQ[row + col]];
	{BZA, BZB}
];*)

(* Return the energy and weight lists for computing local G.F. *)
GetLatticeEnergies[HalfBandwidths_, \[Delta]_, LatticeType_, LatticeDim_, NumberOfPoints_] := Module[
	{energies, weights, Norb = Length[HalfBandwidths], d\[Epsilon], LE, BZ},
	Which[
		LatticeType == "Bethe" && LatticeDim == Infinity,
		(* initialize energies list of rank LE x Norb x Norb *)
		energies = ConstantArray[0, {NumberOfPoints, Norb, Norb}];
		d\[Epsilon] = 2.0*HalfBandwidths[[1]]/(NumberOfPoints-1);
		(* energies are sampled at regular intervals *)
		energies[[All, 1, 1]] = Table[\[Epsilon], {\[Epsilon], -HalfBandwidths[[1]], HalfBandwidths[[1]], d\[Epsilon]}];
		(* weights correspond to the DoS of orb=1 *)
		weights = d\[Epsilon] * DoSBethe[energies[[All,1,1]], HalfBandwidths[[1]]];
		(* fill up orbital diagonal energies starting from the last to the first orbital (energies[[All, 1,1]] must be updated at the very end) *)
		Do[
			energies[[All, orb, orb]] = (HalfBandwidths[[orb]]/HalfBandwidths[[1]]) * energies[[All, 1, 1]] + \[Delta][[orb]];
		, {orb, Norb, 1, -1}];,
	(* ---------------------------------------------------------------------- *)
		LatticeType == "Hypercubic" && LatticeDim != Infinity,
		(* number of points per lattice direction *)
		LE = Floor[NumberOfPoints^(1./LatticeDim)];
		(* get the Brillouin Zone *)
		BZ = BrillouinZone[LE, LatticeDim, LatticeType];
		(* initialize energies list of rank LE x Norb x Norb *)
		energies = ConstantArray[0, {Length[BZ], Norb, Norb}];
		(* equal weights to all the energies since we are sampling the Brillouin zone *)
		weights = ConstantArray[1./(Length[BZ]), Length[BZ]];
		Do[
			energies[[All, orb, orb]] = DispersionHypercubic[BZ, HalfBandwidths[[orb]]] + \[Delta][[orb]];
		, {orb, Norb}];,
	(* ---------------------------------------------------------------------- *)
		LatticeType != "Hypercubic",
		Print["WARNING!!! If e(-k) != e(k), the local Green function is WRONG in case Nambu formalism is needed."]
	];
	{energies, weights}
];

(* Return the energy and weight lists for computing local G.F. when EdMode = "Raman" and spin states are not symmetric *)
GetLatticeEnergiesRaman[HalfBandwidths_, \[Delta]_, M_, \[Gamma]_, LatticeType_, LatticeDim_, NumberOfPoints_] := Module[
	{LE, BZ, energies, weights, f = Length[M[[1]]], Norb = Length[HalfBandwidths], P, Pdg, eigvecs},
	Which[
		LatticeType == "Hypercubic" && LatticeDim != Infinity,
		(* number of points per lattice direction *)
		LE = Floor[NumberOfPoints^(1./LatticeDim)];
		(* get the Brillouin Zone *)
		BZ = BrillouinZone[LE, LatticeDim, LatticeType];
		(* initialize energies list of rank LE x Norb x Norb x f x f *)
		energies = ConstantArray[ConstantArray[0.0, {f, f}], {Length[BZ], Norb, Norb}];
		(* equal weights to all the energies since we are sampling the Brillouin zone *)
		weights = ConstantArray[1./(Length[BZ]), Length[BZ]];
		Do[
			(* get the unitary matrix that diagonalizes M: P.M.Pdg = \[CapitalLambda], where \[CapitalLambda] is diagonal *)
			(*eigvecs = Last[ SortBy[Eigensystem[M[[orb]]]\[Transpose], First]\[Transpose] ]; (* list of eigenvectors sorted by eigenvalues (from lower to higher) *)
			Pdg = (Normalize[#] &/@ eigvecs)\[Transpose];
			P = ConjugateTranspose[Pdg];*)
			energies[[All, orb, orb]] = 
				(DispersionHypercubicRaman[#, HalfBandwidths[[orb]], M[[orb]], \[Gamma][[orb]]] + \[Delta][[orb]]*IdentityMatrix[f]) &/@ BZ
		, {orb, Norb}];
	];
	{energies, weights}
];

(* Momentum resolved lattice Hamiltonian H_k in presence of 2 sublattices for Raman mode *)
GetLatticeEnergiesRamanSublattices[HalfBandwidths_, \[Delta]_, M_, \[Gamma]_, LatticeType_, LatticeDim_, NumberOfPoints_] := Module[
	{LE, MBZ, energies, weights, f = Length[M[[1]]], Norb = Length[HalfBandwidths], W = HalfBandwidths},
	If[LatticeType != "Hypercubic", Print["Lattice type not supported."]; ];
	(* number of points per lattice direction *)
	LE = Floor[NumberOfPoints^(1./LatticeDim)];
	(* get the Brillouin Zone *)
	MBZ = MagneticBrillouinZone[LE, LatticeDim, LatticeType];
	(* initialize energies list of rank LE x Norb x Norb x f x f *)
	energies = ConstantArray[ConstantArray[0.0, {2f, 2f}], {Length[MBZ], Norb, Norb}];
	(* equal weights to all the energies since we are sampling the Brillouin zone *)
	weights = ConstantArray[1./(Length[MBZ]), Length[MBZ]];
	(* fill up matrix elements *)
	Which[
		LatticeDim == 1,
		Do[
			energies[[All, orb, orb]] = (ArrayFlatten[{
				{
					M[[orb]] + \[Delta][[orb]]*IdentityMatrix[f], 
					-W[[orb]]*DiagonalMatrix[
						Table[
							(Exp[I*m*\[Gamma][[orb,1]]] + Exp[-2I*#[[1]]]*Exp[-I*m*\[Gamma][[orb,1]]])
						, {m, (f-1)/2, -(f-1)/2, -1}]
					]
				},
				{
					-W[[orb]]*DiagonalMatrix[
						Table[
							(Exp[-I*m*\[Gamma][[orb,1]]] + Exp[2I*#[[1]]]*Exp[I*m*\[Gamma][[orb,1]]])
						, {m, (f-1)/2, -(f-1)/2, -1}]
					], 
					M[[orb]] + \[Delta][[orb]]*IdentityMatrix[f]
				}
			}]) &/@ MBZ
		, {orb, Norb}];,
	(* -------------------------------------------- *)
		LatticeDim == 2,
		Do[
			energies[[All, orb, orb]] = (ArrayFlatten[{
				{
					(* A-A block *)
					M[[orb]] + \[Delta][[orb]]*IdentityMatrix[f], 
					(* A-B block *)
					-W[[orb]]*DiagonalMatrix[Table[
						Exp[I*m*\[Gamma][[orb,1]]] 
						+ Exp[-I*m*\[Gamma][[orb,1]]]*Exp[-2I*#[[1]]] 
						+ Exp[I*m*\[Gamma][[orb,2]]]*Exp[I*(-#[[1]]+#[[2]])] 
						+ Exp[-I*m*\[Gamma][[orb,2]]]*Exp[-I*(#[[1]]+#[[2]])]
					, {m, (f-1)/2, -(f-1)/2, -1}]]
				},
				{
					(* B-A block *)
					-W[[orb]]*DiagonalMatrix[Table[
						Exp[-I*m*\[Gamma][[orb,1]]] 
						+ Exp[I*m*\[Gamma][[orb,1]]]*Exp[2I*#[[1]]] 
						+ Exp[-I*m*\[Gamma][[orb,2]]]*Exp[-I*(-#[[1]]+#[[2]])] 
						+ Exp[I*m*\[Gamma][[orb,2]]]*Exp[I*(#[[1]]+#[[2]])]
					, {m, (f-1)/2, -(f-1)/2, -1}]], 
					(* B-B block *)
					M[[orb]] + \[Delta][[orb]]*IdentityMatrix[f]
				}
			}]) &/@ MBZ
		, {orb, Norb}];
	];
	{energies, weights}
];


(*           LOCAL GREEN FUNCTION         *)
(* when EdMode = "Normal" *)
LocalGreenFunctionNormal = Compile[{
	{Energies,_Real,1}, {weights,_Real,1}, {\[Mu], _Real}, {\[CapitalSigma], _Complex, 1}, {zlist, _Complex, 1}
	},
	With[
		{LE = Length[Energies]},
		(* compute Gloc *)
		Sum[
			weights[[i]]/(zlist + \[Mu] - Energies[[i]] - \[CapitalSigma])
		, {i, LE}]
	],
	RuntimeAttributes->{Listable}, Parallelization->True
];

(* when EdMode == "Superc". This is based on analytic inversion of the matrix *)
LocalGreenFunctionSupercAnalytic = Compile[{
	{Energies,_Real,1}, {weights, _Real,1},{\[Mu], _Real}, {\[CapitalSigma], _Complex, 3}, {zlist, _Complex, 1}
	},
	Module[
		{LE = Length[Energies], LocalGF = 0*\[CapitalSigma]},
		(* compute normal component *)
		LocalGF[[All,1,1]] = Sum[
			weights[[i]]*(- zlist + \[Mu] - Conjugate@\[CapitalSigma][[All,1,1]] - Energies[[i]])/(Abs[zlist + \[Mu] - \[CapitalSigma][[All,1,1]] - Energies[[i]]]^2 + Abs[\[CapitalSigma][[All,1,2]]]^2)
		,{i, LE}];
		(* compute anomalous component *)
		LocalGF[[All,1,2]] = -\[CapitalSigma][[All,1,2]]*Sum[
			weights[[i]]/(Abs[zlist + \[Mu] - \[CapitalSigma][[All,1,1]] - Energies[[i]]]^2 + Abs[\[CapitalSigma][[All,1,2]]]^2)
		,{i, LE}];
		(* compute the full matrix *)
		LocalGF[[All,2,2]] = - Conjugate[LocalGF[[All,1,1]]];
		LocalGF[[All,2,1]] = Conjugate[LocalGF[[All,1,2]]];
		LocalGF
	],
	RuntimeAttributes->{Listable}, Parallelization->True
];
(* use explicit inversion of 2x2 matrices *)
LocalGreenFunctionSuperc = Compile[{
	{Energies,_Real,1}, {weights, _Real,1}, {\[Mu], _Real}, {\[CapitalSigma], _Complex, 3}, {zlist, _Complex, 1}
	},
	Module[
		{LE = Length[Energies], NMatsubara = Length[zlist]},
		Total @ TwoByTwoInverse[
				Table[(1./weights[[i]]) * (
					(IdentityMatrix[2] * #) &/@ zlist +
					ConstantArray[
						(\[Mu] - Energies[[i]]) * DiagonalMatrix[{1., -1.}]
					, NMatsubara] - \[CapitalSigma])
				, {i, 2, LE-1}]
			]
	],
	RuntimeAttributes->{Listable}, Parallelization->True
];

(* when EdMode == "Raman" *)
LocalGreenFunctionRaman = Compile[{
	{Energies,_Real,3}, {weights, _Real,1}, {\[Mu], _Real}, {\[CapitalSigma], _Complex, 3}, {zlist, _Complex, 1}
	},
	Module[
		{LE = Length[Energies], NMatsubara = Length[zlist], f = Length[Energies[[1]]]},
		If[f == 2, (* in this case use faster inversion algorithm *)
			Total @ TwoByTwoInverse[
				Table[(1./weights[[i]]) * (
					(IdentityMatrix[f] * #) &/@ zlist +
					ConstantArray[
						(\[Mu]*IdentityMatrix[f] - Energies[[i]])
					, NMatsubara] - \[CapitalSigma])
				, {i, 1, LE}]
			],
		(* else, if f>2 use built in inversion procedure *)
			Total @ Table[
				weights[[i]] * Inverse[#]&/@(
					(IdentityMatrix[f]*#)&/@zlist +
					ConstantArray[
						\[Mu] * IdentityMatrix[f] - Energies[[i]]
					, NMatsubara] - \[CapitalSigma]
				)
			,{i,LE}]
		]
	],
	RuntimeAttributes->{Listable}, Parallelization->True
];

(* when EdMode == "Raman" and there is sublattice structure *)
LocalGreenFunctionRamanSublattices = Compile[{
	{Energies,_Complex,3}, {weights, _Real,1}, {\[Mu], _Real}, {h,_Real}, {\[CapitalSigma], _Complex, 3}, {zlist, _Complex, 1}
	},
	Module[
		{LE = Length[Energies], NMatsubara = Length[zlist], dim = Length[Energies[[1]]], Gloc},
		(* now use this matrix to compute Gloc *)
		Gloc = Total @ Table[
			weights[[i]] * Inverse[#]&/@(
				(IdentityMatrix[dim]*#)&/@zlist +
				ConstantArray[
					\[Mu] * IdentityMatrix[dim] - Energies[[i]] - DiagonalMatrix[{h,-h,-h,h}]
				, NMatsubara] - \[CapitalSigma]
			)
		, {i,LE}];
		(* return only block diagonal part, i.e. {Gloc_sublatticeA, Gloc_sublatticeB} *)
		{Gloc[[All, 1;;2, 1;;2]], Gloc[[All, 3;;4, 3;;4]]}
	],
	RuntimeAttributes->{Listable}, Parallelization->True
];


(* when EdMode = "InterorbNormal" *)
LocalGreenFunctionInterorbNormal = Compile[{
	{Energies,_Real,3}, {weights,_Real,1}, {\[Mu], _Real}, {\[CapitalSigma], _Complex, 3}, {zlist, _Complex, 1}
	},
	With[
		{LE = Length[Energies], Norb = Length[Energies[[1]]], NMatsubara = Length[zlist]},
		(* compute Gloc *)	
		Total@Table[
			weights[[i]]*Inverse[#]&/@(
				(IdentityMatrix[Norb]*#)&/@(zlist+\[Mu]) - ConstantArray[Energies[[i]],NMatsubara] - \[CapitalSigma]
			)
		,{i,LE}]
	],
	RuntimeAttributes->{Listable}, Parallelization->True
];

(* when EdMode = "FullSuperc" *)
LocalGreenFunctionFullSuperc = Compile[{
	{Energies,_Real,3}, {weights,_Real,1}, {\[Mu], _Real}, {\[CapitalSigma], _Complex, 3}, {zlist, _Complex, 1}
	},
	With[
		{LE = Length[Energies], Norb = Round[Length[Energies[[1]]]/2], NMatsubara = Length[zlist]},
		(* compute Gloc *)	
		Total@Table[
			weights[[i]] * Inverse[#]&/@(
				(IdentityMatrix[2*Norb]*#)&/@zlist +
				ConstantArray[
					\[Mu] * DiagonalMatrix[Table[(-1)^(j+1),{j,2*Norb}]] - Energies[[i]]
				, NMatsubara] - \[CapitalSigma]
			)
		,{i,LE}]
	],
	RuntimeAttributes->{Listable}, Parallelization->True
];

(* Local Green function sbatch *)
LocalGreenFunction[LatticeEnergies_, weights_, \[Mu]_, \[CapitalSigma]_, zlist_, EdMode_, SublatticesQ_, OptionsPattern[]] := Module[
	{energies, LE = Length[LatticeEnergies], Norb = Length[LatticeEnergies[[1]]], \[CapitalSigma]matrix},
	Which[
		EdMode == "Normal" && !SublatticesQ,
		(* the input will be {{e1},{e2},...} but we want {e1, e2, e3, ...} *)
		LocalGreenFunctionNormal[Flatten[LatticeEnergies], weights, \[Mu], \[CapitalSigma], zlist],
	(* ------------------------------------------------------------------- *)
		EdMode == "Superc" && !SublatticesQ,
		If[Re[zlist[[-1]]] == 0.0, (* if we are using Matsubara frequencies, use a shortcut *)
			LocalGreenFunctionSupercAnalytic[Flatten[LatticeEnergies], weights, \[Mu], \[CapitalSigma], zlist],
		(* else if we are using real frequencies let's be safe and invert the matrix *)
			LocalGreenFunctionSuperc[Flatten[LatticeEnergies], weights, \[Mu], \[CapitalSigma], zlist]
		],
	(* ------------------------------------------------------------------- *)
		EdMode == "Raman" && !SublatticesQ,
		(* the input LatticeEnergies will be a list of fxf matrices (orbital indexes will be specified) *)
		LocalGreenFunctionRaman[LatticeEnergies, weights, \[Mu], \[CapitalSigma], zlist],
	(* ------------------------------------------------------------------- *)
		EdMode == "Raman" && SublatticesQ,
		(* the input LatticeEnergies will be a list of 2fx2f matrices (orbital indexes will be specified) *)
		LocalGreenFunctionRamanSublattices[LatticeEnergies, weights, \[Mu], OptionValue["StaggeredMagneticField"], \[CapitalSigma], zlist],
	(* ------------------------------------------------------------------- *)
		EdMode == "InterorbNormal",
		(* in this situation the input tensor has the correct shape *)
		LocalGreenFunctionInterorbNormal[LatticeEnergies, weights, \[Mu], \[CapitalSigma], zlist],
	(* ------------------------------------------------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		energies = ConstantArray[0, {LE, 2*Norb, 2*Norb}];
		(* diagonal part of energy tensor in Nambu representation *)
		Do[
			(* top left element of diagonal blocks *)
			energies[[All, 2(orb-1)+1, 2(orb-1)+1]] = LatticeEnergies[[All, orb, orb]];
			(* bottom right element of diagonal blocks *)
			energies[[All, 2(orb-1)+2, 2(orb-1)+2]] = - LatticeEnergies[[All, orb, orb]];
		, {orb, Norb}];
		LocalGreenFunctionFullSuperc[energies, weights, \[Mu], \[CapitalSigma], zlist]
	]
];
Options[LocalGreenFunction] = {"StaggeredMagneticField" -> 0.0}


End[]

EndPackage[]
