(* ::Package:: *)

BeginPackage["Lattices`"]


DoSBethe::usage = "."

BrillouinZone::usage = "."

GetLatticeEnergies::usage = "."

DispersionHypercubic::usage = "."

LocalGreenFunction::usage = "LocalGreenFunction[LatticeEnergies_, weights_, \[Mu]_, \[CapitalSigma]_, zlist_, EdMode_]"


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

(* returns a list of k points in the 1st BZ *)
BrillouinZone[LE_, d_, OptionsPattern[]] := With[
	{dk = 2.Pi/LE, Lattice = OptionValue[Lattice]},
	Which[
		Lattice == "Hypercubic",
		Flatten[
			Outer[{##}&,##]&@@
				ConstantArray[
					Table[k, {k, -1.*Pi, 1.*Pi-dk, dk}],
					d
				],
		d-1]
	]
];
Options[BrillouinZone] = {Lattice -> "Hypercubic"};

(* Return the energy and weight lists for computing local G.F. *)
GetLatticeEnergies[HalfBandwidths_, LatticeType_, LatticeDim_, NumberOfPoints_] := Module[
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
		(* fill up orbital diagonal energies *)
		Do[
			energies[[All, orb, orb]] = (HalfBandwidths[[orb]]/HalfBandwidths[[1]]) * energies[[All, 1, 1]];
		, {orb, 2, Norb}];,
	(* ---------------------------------------------------------------------- *)
		LatticeType == "Hypercubic" && LatticeDim != Infinity,
		(* number of points per lattice direction *)
		LE = Floor[NumberOfPoints^(1./LatticeDim)];
		(* get the Brillouin Zone *)
		BZ = BrillouinZone[LE, LatticeDim, Lattice -> "Hypercubic"];
		(* initialize energies list of rank LE x Norb x Norb *)
		energies = ConstantArray[0, {Length[BZ], Norb, Norb}];
		(* equal weights to all the energies since we are sampling the Brillouin zone *)
		weights = ConstantArray[1./(Length[BZ]), Length[BZ]];
		Do[
			energies[[All, orb, orb]] = DispersionHypercubic[BZ, HalfBandwidths[[orb]]];
		, {orb, Norb}];,
	(* ---------------------------------------------------------------------- *)
		LatticeType != "Hypercubic",
		Print["WARNING!!! If e(-k) != e(k), the local Green function is WRONG in case Nambu formalism is needed."]
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

(* when EdMode == "Superc" *)
LocalGreenFunctionSuperc = Compile[{
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
LocalGreenFunction[LatticeEnergies_, weights_, \[Mu]_, \[CapitalSigma]_, zlist_, EdMode_] := Module[
	{energies, LE = Length[LatticeEnergies], Norb = Length[LatticeEnergies[[1]]]},
	Which[
		EdMode == "Normal",
		energies = Flatten[LatticeEnergies]; (* the input will be {{e1},{e2},...} but we want {e1, e2, e3, ...} *)
		LocalGreenFunctionNormal[energies, weights, \[Mu], \[CapitalSigma], zlist],
	(* ------------------------------------------------------------------- *)
		EdMode == "Superc",
		(* --- this code is ok for computing LocalG via inverse matrix
		energies = ConstantArray[0, {LE, 2, 2}]; (* the input will be {{e1},{e2},...} but we want the 2x2 Nambu representation *)
		energies[[All, 1, 1]] = Flatten[LatticeEnergies];
		energies[[All, 2, 2]] = - energies[[All, 1, 1]];
		*)
		energies = Flatten[LatticeEnergies];
		LocalGreenFunctionSuperc[energies, weights, \[Mu], \[CapitalSigma], zlist],
	(* ------------------------------------------------------------------- *)
		EdMode == "InterorbNormal",
		energies = LatticeEnergies; (* in this situation the input tensor has the correct shape *)
		LocalGreenFunctionInterorbNormal[energies, weights, \[Mu], \[CapitalSigma], zlist],
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



End[]

EndPackage[]
