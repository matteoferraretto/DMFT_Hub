(* ::Package:: *)

$Path = Join[$Path, {"C:\\Users\\matte\\Desktop\\Mathematica_package\\"}];
(* Import the code *)
<<"DMFT.wl";


GetLatticeEnergiesRaman[HalfBandwidths_, \[Delta]_, M_, \[Gamma]_, u_, LatticeType_, LatticeDim_, NumberOfPoints_] := Module[
	{LE, BZ, energies, weights, f = Length[M], Norb = Length[HalfBandwidths], P, Pdg},
	Which[
		LatticeType == "Hypercubic" && LatticeDim != Infinity,
		(* number of points per lattice direction *)
		LE = Floor[NumberOfPoints^(1./LatticeDim)];
		(* get the Brillouin Zone *)
		BZ = BrillouinZone[LE, LatticeDim, Lattice -> "Hypercubic"];
		(* get the unitary matrix that diagonalizes M: P.M.Pdg = \[CapitalLambda], where \[CapitalLambda] is diagonal *)
		Pdg = Normalize[#] &/@ Eigenvectors[M]\[Transpose];
		P = ConjugateTranspose[Pdg];
		(* initialize energies list of rank LE x Norb x Norb x f x f *)
		energies = ConstantArray[ConstantArray[0, {f, f}], {Length[BZ], Norb, Norb}];
		(* equal weights to all the energies since we are sampling the Brillouin zone *)
		weights = ConstantArray[1./(Length[BZ]), Length[BZ]];
		Do[
			energies[[All, orb, orb]] = (P . # . Pdg) &/@ ((DispersionHypercubicRaman[#, HalfBandwidths[[orb]], f, \[Gamma], u] + \[Delta][[orb]]*IdentityMatrix[f] + M) &/@ BZ);
		, {orb, Norb}];
	];
	{energies, weights}
];


L = 3; Norb = 2; f = 2;
LatticeType = "Hypercubic"; LatticeDim = 1; NumberOfPoints = 400;
W = {1., 1.};
\[Delta] = {1., 0.};
M = 0.5 * PauliMatrix[1];
h = Eigenvalues[M];
\[Gamma] = Pi/4.;
u = {1., 0.};
BZ = BrillouinZone[NumberOfPoints, LatticeDim];
{LatticeEnergies, LatticeWeights} = GetLatticeEnergiesRaman[W, \[Delta], M, \[Gamma], u, LatticeType, LatticeDim, NumberOfPoints];


Manipulate[

GraphicsColumn[{
	Show[
		ListPlot[
		{Flatten[BZ], LatticeEnergies[[All, orb, orb, \[Sigma], \[Sigma]]]}\[Transpose]
		],
		Plot[-2.*Cos[\[Gamma]/2]*Cos[k] + \[Delta][[orb]] + h[[\[Sigma]]], {k, -Pi, Pi}, PlotStyle->Directive[Orange, Dashing[.05]]]
	],

	Show[
		ListPlot[
		{Flatten[BZ], LatticeEnergies[[All, orb, orb, \[Sigma], Mod[\[Sigma],f]+1 ]]}\[Transpose]
		],
		Plot[2.*Sin[\[Gamma]/2]*Sin[k], {k, -Pi, Pi}, PlotStyle->Directive[Orange, Dashing[.05]]]
	]
}]

, {orb, Range[Norb]}, {\[Sigma], Range[f]}]


Pdg = Normalize[#] &/@ Eigenvectors[PauliMatrix[1]]\[Transpose];
P = ConjugateTranspose[Pdg];

FullSimplify[(P . ({
{Cos[k - g/2], 0},
{0, Cos[k + g/2]} 
}+ a*PauliMatrix[1]) . Pdg)]//MatrixForm
