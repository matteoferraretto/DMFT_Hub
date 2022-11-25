(* ::Package:: *)

BeginPackage["Observables`", {"Sectors`", "Hamiltonian`", "ImpurityGreenFunction`", "Lattices`"}]


Density::usage = "Density[L, f, Norb, j, \[Sigma], orb, Sectors, EgsSectorList, GsSectorList, T]"

SquareDensity::usage = "SquareDensity[L, f, Norb, {i, j}, {\[Sigma]1, \[Sigma]2}, {orb1, orb2}, Sectors, EgsSectorList, GsSectorList, T]"

CdgCdg::usage = "CdgCdg[L, f, Norb, {i,j}, {\[Sigma]1,\[Sigma]2}, {orb1,orb2}, Sectors, EgsSectorList, GsSectorList, T]"

CdgC::usage = "CdgC[L, f, Norb, {i,j}, {\[Sigma]1,\[Sigma]2}, {orb1,orb2}, Sectors, EgsSectorList, GsSectorList, T]"

QuasiparticleWeight::usage = "QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode] computes the quasiparticle weight z from the Self-energy. "

OrderParameter::usage = "OrderParameter[InverseG, TMats] returns the superconductive order parameter computed from the Green function. "

SuperfluidStiffness::usage = "SuperfluidStiffness[DBethe, \[CapitalSigma], i\[Omega]] returns the superfluid stiffness computed from the Green function. "

KineticEnergy::usage = "KineticEnergy[DBethe, \[Mu], \[CapitalSigma], i\[Omega], EdMode] computes the Kinetic energy (expectation value of non-local impurity Hamiltonian). "

SpectralFunction::usage = "SpectralFunction[L, f, Norb, \[Sigma], orb, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega], \[Eta]]"


Begin["Private`"]

Print["Package Observables` loaded successfully."];

(*                      OPERATORS                      *)
(* Density *)
Density[L_, f_, Norb_, j_, \[Sigma]_, orb_, Sectors_, EgsSectorList_, GsSectorList_, T_, OptionsPattern[]] := Module[
	{Egs, Gs, GsQns, GsSectorIndex, \[Epsilon] = OptionValue[DegeneracyThreshold], num, \[Psi], density},
	If[T == 0,
		Egs = Min[Flatten[EgsSectorList]];(* ground state energy (lowest of all the sectors) *)
		GsSectorIndex = Position[
			EgsSectorList,
			_?((Abs[# - Egs] < \[Epsilon])&)
		]; (* sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy *) 
		density = Sum[
			\[Psi] = Sectors[[index[[1]]]];
			num = n[L, f, Norb, j, \[Sigma], orb, \[Psi]]; (* local density *)
			num . Abs[GsSectorList[[##]]&@@index]^2
		, {index, GsSectorIndex}]
	];
	density/Length[GsSectorIndex]
];
Options[Density] = {DegeneracyThreshold -> 1.*10^(-9)}

(* Density-Density *)
SquareDensity[L_, f_, Norb_, {i_,j_}, {\[Sigma]1_,\[Sigma]2_}, {orb1_,orb2_}, Sectors_, EgsSectorList_, GsSectorList_, T_, OptionsPattern[]] := Module[
	{Egs, Gs, GsQns, GsSectorIndex, \[Epsilon] = OptionValue[DegeneracyThreshold], num, \[Psi], squaredensity},
	If[T == 0,
		Egs = Min[Flatten[EgsSectorList]];(* ground state energy (lowest of all the sectors) *)
		GsSectorIndex = Position[
			EgsSectorList,
			_?((Abs[# - Egs] < \[Epsilon])&)
		]; (* sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy *) 
		squaredensity = Sum[
			\[Psi] = Sectors[[index[[1]]]];
			num = n[L, f, Norb, i, \[Sigma]1, orb1, \[Psi]] * n[L, f, Norb, j, \[Sigma]2, orb2, \[Psi]]; (* n squared *)
			num . Abs[GsSectorList[[##]]&@@index]^2
		, {index, GsSectorIndex}]
	];
	squaredensity/Length[GsSectorIndex]
];
Options[SquareDensity] = {DegeneracyThreshold -> 1.*10^(-9)}

(* < cdg_(i \[Sigma]1 orb1) cdg_(j \[Sigma]2 orb2) > *)
CdgCdg[L_, f_, Norb_, {i_,j_}, {\[Sigma]1_,\[Sigma]2_}, {orb1_,orb2_}, Sectors_, EgsSectorList_, GsSectorList_, T_, OptionsPattern[]] := Module[
	{Egs, Gs, GsQns, GsSectorIndex, \[Epsilon] = OptionValue[DegeneracyThreshold], cdgcdg, dim, gs, \[Psi], \[Psi]1, \[Chi], rows, cols, pos, \[CapitalSigma], dispatch, \[Phi]},
	If[T == 0,
		Egs = Min[Flatten[EgsSectorList]];(* ground state energy (lowest of all the sectors) *)
		GsSectorIndex = Position[
			EgsSectorList,
			_?((Abs[# - Egs] < \[Epsilon])&)
		]; (* sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy *) 
		\[Phi] = Sum[
			gs = GsSectorList[[##]]&@@index;
			\[Psi] = Sectors[[index[[1]]]];
			dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1]];
			dim = Length[\[Psi]];
			\[Psi]1 = CreatePairSelect[L, f, i, j, \[Sigma]1, \[Sigma]2, orb1, orb2, \[Psi]];
			cdgcdg = SparseArray[{}, {dim,dim}];
			If[Length[\[Psi]1] != 0,
				\[Chi] = CreatePair[L, f, i, j, \[Sigma]1, \[Sigma]2, orb1, orb2, \[Psi]1];
				rows = \[Chi]/.dispatch;(* *)cols = \[Psi]1/.dispatch;(* *)pos = {rows,cols}\[Transpose];
				\[CapitalSigma] = CCSign[L, f, {i,j}, {\[Sigma]1,\[Sigma]2}, {orb1,orb2}, \[Psi]1];
				cdgcdg += SparseArray[pos -> \[CapitalSigma], {dim,dim}];
			];
			Conjugate[gs] . (cdgcdg . gs)
		, {index, GsSectorIndex}]
	];
	\[Phi]/Length[GsSectorIndex]
];
Options[CdgCdg] = {DegeneracyThreshold -> 1.*10^(-9)}

(* < cdg_(i \[Sigma]1 orb1) c_(j \[Sigma]2 orb2) > *)
CdgC[L_, f_, Norb_, {i_,j_}, {\[Sigma]1_,\[Sigma]2_}, {orb1_,orb2_}, Sectors_, EgsSectorList_, GsSectorList_, T_, OptionsPattern[]] := Module[
	{Egs, Gs, GsQns, GsSectorIndex, \[Epsilon] = OptionValue[DegeneracyThreshold], cdgc, dim, gs, \[Psi], \[Psi]1, \[Chi], rows, cols, pos, \[CapitalSigma], dispatch, hop},
	If[T == 0,
		Egs = Min[Flatten[EgsSectorList]];(* ground state energy (lowest of all the sectors) *)
		GsSectorIndex = Position[
			EgsSectorList,
			_?((Abs[# - Egs] < \[Epsilon])&)
		]; (* sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy *) 
		hop = Sum[
			gs = GsSectorList[[##]]&@@index;
			\[Psi] = Sectors[[index[[1]]]];
			dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1]];
			dim = Length[\[Psi]];
			\[Psi]1 = HopSelect[L, f, i, j, \[Sigma]1, \[Sigma]2, orb1, orb2, \[Psi]];
			cdgc = SparseArray[{}, {dim,dim}];
			If[Length[\[Psi]1] != 0,
				\[Chi] = Hop[L, f, i, j, \[Sigma]1, \[Sigma]2, orb1, orb2, \[Psi]1];
				rows = \[Chi]/.dispatch;(* *)cols = \[Psi]1/.dispatch;(* *)pos = {rows,cols}\[Transpose];
				\[CapitalSigma] = CCSign[L, f, {i,j}, {\[Sigma]1,\[Sigma]2}, {orb1,orb2}, \[Psi]1];
				If[j < i, \[CapitalSigma] = -\[CapitalSigma]];(* if j=L, we are applying cdg_L c_1. When moving cdg_L before position L, it jumps over c_1 and the sign changes! *)
				cdgc += SparseArray[pos -> \[CapitalSigma], {dim,dim}];
			];
			Conjugate[gs] . (cdgc . gs)
		, {index, GsSectorIndex}]
	];
	hop/Length[GsSectorIndex]
];
Options[CdgC] = {DegeneracyThreshold -> 1.*10^(-9)}

(* Spectral function *)
SpectralFunction[L_, f_, Norb_, \[Sigma]_, orb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, \[Omega]_, \[Eta]_] := Module[
	{spectralfunction},
	Which[
		EdMode == "Normal",
		spectralfunction = -(1./Pi) * Im[Mean[MapApply[
			GreenFunctionImpurity[L, f, Norb, \[Sigma], orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega]+I*\[Eta]]&,
			{Gs, GsQns}\[Transpose]
		]]],
(* ---------------------------------------------- *)
		EdMode == "Superc",
		spectralfunction = -(1./Pi) * (1./f) * Im[ Tr[#] &/@ 
			Mean[MapApply[
				GreenFunctionImpurityNambu[L, f, Norb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega]+I*\[Eta], Orb -> orb]&,
				{Gs, GsQns}\[Transpose]
			]]]
	];
	{\[Omega], spectralfunction}\[Transpose]
];

(* Quasiparticle weight *)
\[NonBreakingSpace]QuasiparticleWeight[\[CapitalSigma]_, i\[Omega]_, EdMode_, OptionsPattern[]] := Module[
	{Selfenergy, data, a, z, cutoff = OptionValue[FitCutoff], orb = OptionValue[Orb]},
	Selfenergy = Which[
		EdMode == "Normal", \[CapitalSigma], 
		EdMode == "Superc", \[CapitalSigma][[All,1,1]],
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc", \[CapitalSigma][[All, 2(orb-1)+1, 2(orb-1)+1]]
	];
	data = ({Im[i\[Omega]], Im[Selfenergy]}\[Transpose])[[;;cutoff]];
	a = Fit[data, {x}, x]/x;
	z = 1./(1.-a)
];
Options[QuasiparticleWeight] = {FitCutoff -> 50, Orb -> 1};

(* \[Phi]: superconductive order parameter computed through the Green function *)
OrderParameter[InverseG_, TMats_] := With[
	{G = Inverse[#] &/@ InverseG},
	(* the factor 2 comes from the fact that we have to account for negative Matsubara frequencies, the Re[] is to suppress a tiny imaginary part *)
	- 2.0 * TMats * Re @ Total[G[[All, 1, 2]]] 
];



(* ------------ in progress ----------- *)
(* Superfluid Stiffness *)
SuperfluidStiffness[DBethe_, \[CapitalSigma]_, i\[Omega]_, OptionsPattern[]] := Module[
	{LE = OptionValue[NumberOfPoints], Lattice = OptionValue[Lattice], dim = OptionValue[LatticeDimension], d\[Epsilon], TMats, Ds},
	(* initialize parameters *)
	TMats = (i\[Omega][[2]] - i\[Omega][[1]])/(2*Pi*I);
	Which[
		Lattice == "Bethe",
		d\[Epsilon] = 2.*DBethe/LE;
		(* compute the stiffness *)
		Ds = 4. * TMats * d\[Epsilon] * Total @ Table[
			DoSBethe[\[Epsilon], DBethe] * (* density of states of the Bethe lattice *)
			((DBethe^2-\[Epsilon]^2)/3.) * (* current vertex function for the Bethe lattice *)
			Total[Abs[-\[CapitalSigma][[All,1,2]]/(Abs[i\[Omega]-\[CapitalSigma][[All,1,1]]-\[Epsilon]]^2+Abs[\[CapitalSigma][[All,1,2]]]^2)]^2] (* \!\(
\*SubscriptBox[\(\[Sum]\), \(i\[Omega]\)]\(\(|\)\(F\((i\[Omega])\)\)
\*SuperscriptBox[\(|\), \(2\)]\)\) *)
		,{\[Epsilon], -DBethe, DBethe, d\[Epsilon]}];
	];
Re[Ds]
];
Options[SuperfluidStiffness] = {Lattice -> "Bethe", LatticeDimension -> 2, NumberOfPoints -> 1000};

(* Kinetic energy, i.e. < Subscript[H, non interacting] > *)
KineticEnergy[DBethe_, \[Mu]_, \[CapitalSigma]_, i\[Omega]_, EdMode_] := Module[
	{LE = OptionValue[NumberOfPoints], Lattice = OptionValue[Lattice], dim = OptionValue[LatticeDimension],TMats = (i\[Omega][[2]]-i\[Omega][[1]])/(2*Pi*I), d\[Epsilon], Glattice, \[CapitalSigma]0, Ekin = 0},
	Which[
		EdMode == "Normal", 
		\[CapitalSigma]0 = Last@\[CapitalSigma]; (* Self energy at the last Matsubara frequency *)
		Glattice[\[Epsilon]_] := 1./(i\[Omega] + \[Mu] - \[Epsilon] - \[CapitalSigma]);, (* compute lattice Green function *)
	(* ----------------------------- *)
		EdMode == "Superc", 
		\[CapitalSigma]0 = \[CapitalSigma][[-1,1,1]];
		Glattice[\[Epsilon]_] := (-i\[Omega]+\[Mu]-Conjugate@\[CapitalSigma][[All,1,1]]-\[Epsilon])/(Abs[i\[Omega]+\[Mu]-\[CapitalSigma][[All,1,1]]-\[Epsilon]]^2+Abs[\[CapitalSigma][[All,1,2]]]^2);
	];
	Which[
		Lattice == "Bethe",
		d\[Epsilon] = 2.*DBethe/LE;		
		(* Non vanishing terms of the first row of Eq. 4.12 - KineticEnergy.PDF *)
		Ekin += 4. * TMats * d\[Epsilon] * Total@Table[
			DoSBethe[\[Epsilon], DBethe] * \[Epsilon] *
			Total[Re@Glattice[\[Epsilon]] - (\[Epsilon] + \[CapitalSigma]0)/(i\[Omega]^2)]
		, {\[Epsilon], -DBethe, DBethe, d\[Epsilon]}];(*  Notice that \!\(
\*SubscriptBox[\(\[Sum]\), \(n\)]\(G\((k, 
\*SubscriptBox[\(i\[Omega]\), \(n\)])\)\)\)=2\!\(
\*SubscriptBox[\(\[Sum]\), \(n \[GreaterEqual] 0\)]\(Re[G\((k, 
\*SubscriptBox[\(i\[Omega]\), \(n\)])\)]\)\)  *)
		(* Non vanishing terms of the second row of Eq. 4.12 - KineticEnergy.PDF *)
		Ekin += (1./(2.*TMats)) * d\[Epsilon] * Total@Table[
			DoS[\[Epsilon]]*\[Epsilon]*(-\[Epsilon] - \[CapitalSigma]0)
		,{\[Epsilon], -DBethe, DBethe, d\[Epsilon]}];,
	(* -------------------------------------------------------- *)
		Lattice == "Hypercubic",
		Return[0]
	];
	Re[Ekin]
];
Options[KineticEnergy] = {Lattice -> "Bethe", LatticeDimension -> 2, NumberOfPoints -> 1000};

End[]

EndPackage[]
