(* ::Package:: *)

BeginPackage["Observables`", {"MyLinearAlgebra`", "Sectors`", "Hamiltonian`", "ImpurityGreenFunction`", "Lattices`"}]


Density::usage = "Density[L, f, Norb, j, \[Sigma], orb, Sectors, EgsSectorList, GsSectorList, T]"

SquareDensity::usage = "SquareDensity[L, f, Norb, {i, j}, {\[Sigma]1, \[Sigma]2}, {orb1, orb2}, Sectors, EgsSectorList, GsSectorList, T]"

MomentumDistributedDensityRaman::usage = "MomentumDistributedDensityRaman[kindexes_, Energies_, \[Mu]_, \[CapitalSigma]_, i\[Omega]_]"

FlavorCurrent::usage = "FlavorCurrent[t, \[Gamma], \[Sigma], a, flavordistribution, LatticeType, LatticeDim, NumberOfPoints] gives the a-th spatial component of the flavor current 
associated to the flavor \[Sigma]=1,2,...,f. This function requires as an input ''flavordistribution'', i.e. the flavor-resolved momentum-distributed operator <cdg_k\[Alpha] c_k\[Beta]>, 
which is a list of fxf matrices, each one associated to a given momentum. "

ChiralCurrentNonInteracting::usage = "ChiralCurrentNonInteracting[t_, \[CapitalOmega]_, \[Gamma]_] returns the non interacting value of the chiral current in the thermodynamic limit by using an analytic
expression."

PsidgPsi::usage = "PsidgPsi[kindexes, H, \[Mu], \[CapitalSigma], i\[Omega]] computes < \[Psi]dag_k \[Psi]_k > from the converged self-energy. Should work for every EdMode. "

CdgCdg::usage = "CdgCdg[L, f, Norb, {i,j}, {\[Sigma]1,\[Sigma]2}, {orb1,orb2}, Sectors, EgsSectorList, GsSectorList, T]"

CdgC::usage = "CdgC[L, f, Norb, {i,j}, {\[Sigma]1,\[Sigma]2}, {orb1,orb2}, Sectors, EgsSectorList, GsSectorList, T]"

QuasiparticleWeight::usage = "QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode] computes the quasiparticle weight z from the Self-energy. "

OrderParameter::usage = "OrderParameter[InverseG, TMats] returns the superconductive order parameter computed from the Green function. "

SuperfluidStiffness::usage = "SuperfluidStiffness[DBethe, \[CapitalSigma], i\[Omega]] returns the superfluid stiffness computed from the Green function. "

KineticEnergy::usage = "KineticEnergy[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega], EdMode] computes the Kinetic energy (expectation value of non-local impurity Hamiltonian). "

SpectralFunction::usage = "SpectralFunction[LatticeEnergies_, weights_, \[Mu]_, \[CapitalSigma]_, zlist_, EdMode_, SublatticesQ_]"

MomentumResolvedSpectralFunction::usage = "MomentumResolvedSpectralFunction"


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

(* n(k) - to be tested *)
MomentumDistributedDensityNormal[i_, Energies_, \[Mu]_, \[CapitalSigma]_, i\[Omega]_] := Module[
	{TMats = Re[(i\[Omega][[2]]-i\[Omega][[1]])/(2.*Pi*I)], \[CapitalSigma]0 = Last @ \[CapitalSigma]},
	2.0 * ( 2.0 * TMats * Total[
		(* G_numerical - G_tail *)
		Re[ 1./((# + \[Mu] - Energies[[i]]) &/@ i\[Omega] - \[CapitalSigma]) 
		- (((Energies[[i]] + \[CapitalSigma]0)/(#^2)) &/@ i\[Omega]) ]
		] + 1./2. - (1./(4.*TMats)) * (Energies[[i]] + \[CapitalSigma]0) )
];

(* to be tested *)
MomentumDistributedDensitySuperc[i_, Energies_, \[Mu]_, \[CapitalSigma]_, i\[Omega]_] := Module[
	{TMats = Re[(i\[Omega][[2]]-i\[Omega][[1]])/(2.*Pi*I)], \[CapitalSigma]0 = Last @ \[CapitalSigma]},
	2.0 * ( 2.0 * TMats * Total[
		(* G_numerical - G_tail *)
		Re[ TwoByTwoInverse[(#*IdentityMatrix[2] + \[Mu]*PauliMatrix[3] - Energies[[i]]) &/@ i\[Omega] - \[CapitalSigma]] 
		- (((Energies[[i]] + \[CapitalSigma]0)/(#^2)) &/@ i\[Omega]) ]
		] + 1./2. - (1./(4.*TMats)) * (Energies[[i]] + \[CapitalSigma]0) )
];

(* returns a list of matrices < cdg_k\[Alpha] c_k\[Beta] > where \[Alpha],\[Beta] are flavor indexes *)
(* this is also ok when there are sublattices, in that case \[Alpha]=(sub_1,\[Sigma]) \[Beta]=(sub_2,\[Rho]) *)  
MomentumDistributedDensityRaman[kindexes_, Energies_, \[Mu]_, \[CapitalSigma]_, i\[Omega]_] := Module[
	{TMats = Re[(i\[Omega][[2]]-i\[Omega][[1]])/(2.*Pi*I)], f = Length[\[CapitalSigma][[1]]], NFit = Floor[Length[i\[Omega]]/50], \[CapitalSigma]0, flavdist},
	\[CapitalSigma]0 = Mean[ \[CapitalSigma][[;;-NFit]] ];
	(* <cdg_k\[Alpha] c_k\[Beta]> = TMats \!\(
\*SubscriptBox[\(\[Sum]\), \(i\[Omega]\)]\ \(Exp[\(-i\[Omega]\[Eta]\)]\ G\((k, \ i\[Omega])\)_\[Beta]\[Alpha]\)\) *)
	flavdist = Table[
		TMats * Total[(
			(* G_numerical - G_tail (positive Matsubara frequencies) *)
			Transpose[Inverse[#]] &/@ (((#+\[Mu])*IdentityMatrix[f] - Energies[[i]]) &/@ i\[Omega] - \[CapitalSigma]) 
			- ((Transpose[Energies[[i]] + \[CapitalSigma]0]/(#^2)) &/@ i\[Omega]) 
		)]
	, {i, kindexes}];
	(* add G_numerical - G_tail (negative Matsubara frequencies) *)
	flavdist = flavdist + (ConjugateTranspose[#] &/@ flavdist);
	(* add tail contribution *)
	flavdist += Table[
		(1./2.)*IdentityMatrix[f] - (1./(8.*TMats)) * ((Energies[[i]] + \[CapitalSigma]0)\[Transpose] + Conjugate[Energies[[i]] + \[CapitalSigma]0])
	, {i, kindexes}]
];


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
				If[
					f*(orb1-1)+\[Sigma]1 > f*(orb2-1)+\[Sigma]2 || (orb1 == orb2 && \[Sigma]1 == \[Sigma]2 && i > j),
					\[CapitalSigma] = -\[CapitalSigma]
				];
				cdgcdg += SparseArray[pos -> \[CapitalSigma], {dim,dim}];
			];
			Conjugate[gs] . (cdgcdg . gs)
		, {index, GsSectorIndex}]
	];
	\[Phi]/Length[GsSectorIndex]
];
Options[CdgCdg] = {DegeneracyThreshold -> 1.*10^(-9)};

(* < \[Psi]dag_k \[Psi]_k > computed from the converged self-energy *)
PsidgPsi[kindexes_, H_, \[Mu]_, \[CapitalSigma]_, i\[Omega]_] := Module[
	{TMats = Re[(i\[Omega][[2]]-i\[Omega][[1]])/(2.*Pi*I)], f = Length[\[CapitalSigma][[1]]], NFit = Floor[Length[i\[Omega]]/50], \[CapitalSigma]0, \[Psi]dg\[Psi]},
	\[CapitalSigma]0 = Mean[ \[CapitalSigma][[;;-NFit]] ];
	(* <cdg_k\[Alpha] c_k\[Beta]> = TMats \!\(
\*SubscriptBox[\(\[Sum]\), \(i\[Omega]\)]\ \(Exp[\(-i\[Omega]\[Eta]\)]\ G\((k, \ i\[Omega])\)_\[Beta]\[Alpha]\)\) *)
	\[Psi]dg\[Psi] = Table[
		TMats * Total[(
			(* G_numerical - G_tail (positive Matsubara frequencies) *)
			Transpose[Inverse[#]] &/@ (((#+\[Mu])*IdentityMatrix[f] - H[[k]]) &/@ i\[Omega] - \[CapitalSigma]) 
			- ((Transpose[H[[k]] + \[CapitalSigma]0]/(#^2)) &/@ i\[Omega]) 
		)]
	, {k, kindexes}];
	(* add G_numerical - G_tail (negative Matsubara frequencies) *)
	\[Psi]dg\[Psi] = \[Psi]dg\[Psi] + (ConjugateTranspose[#] &/@ \[Psi]dg\[Psi]);
	(* add tail contribution *)
	\[Psi]dg\[Psi] += Table[
		(1./2.)*IdentityMatrix[f] - (1./(8.*TMats)) * ((H[[k]] + \[CapitalSigma]0)\[Transpose] + Conjugate[H[[k]] + \[CapitalSigma]0])
	, {k, kindexes}]
];

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
				If[Index[L, f, Norb, i, \[Sigma]1, orb1] > Index[L, f, Norb, j, \[Sigma]2, orb2], 
					\[CapitalSigma] = -\[CapitalSigma]
				];
				cdgc += SparseArray[pos -> \[CapitalSigma], {dim,dim}];
			];
			Conjugate[gs] . (cdgc . gs)
		, {index, GsSectorIndex}]
	];
	hop/Length[GsSectorIndex]
];
Options[CdgC] = {DegeneracyThreshold -> 1.*10^(-9)};

(* Spectral function of the impurity problem *)
SpectralFunctionImpurity[L_, f_, Norb_, \[Sigma]_, orb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, \[Omega]_, \[Eta]_] := Module[
	{spectralfunction},
	Which[
		EdMode == "Normal",
		spectralfunction = -(1./Pi) * Im[Mean[MapApply[
			GreenFunctionImpurity[L, f, Norb, {orb, orb}, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega]+I*\[Eta]]&,
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

(* spectral function of the lattice problem *)
SpectralFunction[LatticeEnergies_, weights_, \[Mu]_, \[CapitalSigma]_, zlist_, EdMode_, SublatticesQ_] := Module[
	{spectralfunction},
	Which[
		EdMode == "Normal",
		spectralfunction = -(1./Pi) * Im[
			LocalGreenFunction[LatticeEnergies, weights, \[Mu], \[CapitalSigma], zlist, EdMode, SublatticesQ]
		],
	(* ---------------------------------------------- *)
		EdMode == "Superc",
		spectralfunction = -(1./Pi) * Im[
			LocalGreenFunction[LatticeEnergies, weights, \[Mu], \[CapitalSigma], zlist, EdMode, SublatticesQ][[All, 1, 1]]
		],
	(* ---------------------------------------------- *)
		EdMode == "Raman",
		spectralfunction = -(1./Pi) * Im[ Tr[#] &/@
			LocalGreenFunction[LatticeEnergies, weights, \[Mu], \[CapitalSigma], zlist, EdMode, SublatticesQ]
		] / Length[\[CapitalSigma][[1]]],
	(* ---------------------------------------------- *)
		EdMode == "Magnetic",
		spectralfunction = -(1./Pi) * Im[ Total[#] &/@
			LocalGreenFunction[LatticeEnergies, weights, \[Mu], \[CapitalSigma], zlist, EdMode, SublatticesQ]
		] / Length[\[CapitalSigma][[1]]],
	(* ---------------------------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		spectralfunction = -(1./Pi) * Im[ Tr[#] &/@
			LocalGreenFunction[LatticeEnergies, weights, \[Mu], \[CapitalSigma], zlist, EdMode]
		]
	];
	{Re[zlist] + \[Mu], spectralfunction}\[Transpose]
];

(* A(k, \[Omega]) = (\[Omega] + \[Mu] - \[Epsilon]_k - \[CapitalSigma](\[Omega])^-1 *)
MomentumResolvedSpectralFunction[LatticeEnergies_, \[Mu]_, \[CapitalSigma]_, kindexes_, zlist_, EdMode_] := Module[
	{spectralfunction, energies},
	(* extract energies indicated by kindexes *)
	energies = LatticeEnergies[[kindexes]];
	(* *)
	Which[
		EdMode == "Normal",
		spectralfunction = Table[
			- (1./Pi) * Im[1./(zlist[[i]] + \[Mu] - energies[[j]] - \[CapitalSigma][[i]])]
		, {i, Length[zlist]}, {j, Length[kindexes]}],
	(* ---------------------------------------------- *)
		EdMode == "InterorbNormal",
		Print["Not supported"];
		spectralfunction = ConstantArray[0.0, Length[kindexes] * Length[zlist] ],
	(* ---------------------------------------------- *)	
		EdMode == "Superc",
		spectralfunction = Table[ -(1./Pi) * Im[
			TwoByTwoInverse[ zlist[[i]]*IdentityMatrix[2] + (\[Mu] - energies[[j]])*PauliMatrix[3] - \[CapitalSigma][[i]] ][[1,1]]
		], {i, 1, Length[zlist]}, {j, 1, Length[kindexes]}],
	(* ---------------------------------------------- *)
		EdMode == "Raman",
		With[{f = Length[\[CapitalSigma][[1]]], LE = Length[kindexes], NReal = Length[zlist]},	
			spectralfunction = SparseArray[ Chop[ (* convert to a sparse array to save memory *)
				Table[
					- (1./Pi) * Im[
						If[f == 2,
							TwoByTwoInverse[ ((zlist[[i]] + \[Mu])*IdentityMatrix[f] - energies[[j]] - \[CapitalSigma][[i]]) ],
						(* else *)
							Inverse[ ((zlist[[i]] + \[Mu])*IdentityMatrix[f] - energies[[j]] - \[CapitalSigma][[i]]) ] 
						]
					],
				{i, NReal}, {j, LE}]
			, 10^(-2)] ]
		],
	(* ---------------------------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		Print["Not supported"];
		spectralfunction = ConstantArray[0.0, Length[kindexes] * Length[zlist] ];
	];
	spectralfunction
];


(* Quasiparticle weight *)
\[NonBreakingSpace]QuasiparticleWeight[\[CapitalSigma]_, i\[Omega]_, EdMode_, OptionsPattern[]] := Module[
	{Selfenergy, data, a, z, cutoff = OptionValue[FitCutoff], orb = OptionValue[Orb], f = Length[\[CapitalSigma][[1]]]},
	Selfenergy = Which[
		EdMode == "Normal", \[CapitalSigma], 
		EdMode == "Superc", \[CapitalSigma][[All,1,1]],
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc", \[CapitalSigma][[All, 2(orb-1)+1, 2(orb-1)+1]]
	];
	(* if EdMode == "Raman"/"Magnetic" there are f distinct weights, one per effective flavor *)
	If[EdMode == "Raman",
		z = Table[
			data = ({Im[i\[Omega]], Im[\[CapitalSigma][[All, \[Sigma], \[Sigma]]]]}\[Transpose])[[;;cutoff]];
			a = Fit[data, {x, x^2}, x, "BestFitParameters"][[1]];
			1./(1.-a)
		, {\[Sigma], f}];
		Return[z]; (* quit the function *)
	];
	If[EdMode == "Magnetic",
		z = Table[
			data = ({Im[i\[Omega]], Im[\[CapitalSigma][[All, \[Sigma]]]]}\[Transpose])[[;;cutoff]];
			a = Fit[data, {x, x^2}, x, "BestFitParameters"][[1]];
			1./(1.-a)
		, {\[Sigma], f}];
		Return[z]; (* quit the function *)
	];
	(* in all the other cases it is just one number (orbital asymmetries are accounted at higher level in the code) *)
	data = ({Im[i\[Omega]], Im[Selfenergy]}\[Transpose])[[;;cutoff]];
	a = Fit[data, {x, x^2}, x, "BestFitParameters"][[1]];
	z = 1./(1.-a)
];
Options[QuasiparticleWeight] = {FitCutoff -> 10, Orb -> 1};

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

(* Kinetic energy, i.e. < Subscript[H, non interacting] >
KineticEnergyOld[DBethe_, \[Mu]_, \[CapitalSigma]_, i\[Omega]_, EdMode_] := Module[
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
Options[KineticEnergy] = {Lattice -> "Bethe", LatticeDimension -> 2, NumberOfPoints -> 1000};*)


(* EdMode = "Normal" means that the orbitals are not coupled by tunneling terms *)
KineticEnergyNormal[\[Mu]_, LatticeEnergies_, LatticeWeights_, \[CapitalSigma]_, i\[Omega]_, OrbitalSymmetry_] := Module[
	{Norb, \[CapitalSigma]0, TMats, NMats, LE, G, Ekin},
	Norb = Length[\[CapitalSigma]];
	Ekin = ConstantArray[0., Norb]; (* initialize Norb orbital contributions to the kinetic energy *)
	TMats = Im[ i\[Omega][[2]] - i\[Omega][[1]] ] / (2*Pi);
	NMats = Length[i\[Omega]];
	LE = Length[LatticeWeights];
	Do[
		\[CapitalSigma]0 = Last[\[CapitalSigma][[orb]]];
		(* one factor 2 accounts for spin degeneracy, the other one for sum over negative Matsubara frequencies *)
		Ekin[[orb]] = 2. * 2. * TMats * Total @ Table[
			LatticeWeights[[k]] * (LatticeEnergies[[k, orb, orb]] - \[Mu]) * Total @ Table[
				Re[ 1./( (i\[Omega][[n]] + \[Mu]) - LatticeEnergies[[k, orb, orb]] - \[CapitalSigma][[orb, n]] ) ]
				- (LatticeEnergies[[k, orb, orb]] + \[CapitalSigma]0)/(i\[Omega][[n]]^2)
			, {n, NMats}]
		, {k, LE}] + 
		0.5 * Total @ Table[
			LatticeWeights[[k]] * (LatticeEnergies[[k, orb, orb]] - \[Mu])
		, {k, LE}] - 
		(1./(2.*TMats)) * Total @ Table[
			LatticeWeights[[k]] * (LatticeEnergies[[k, orb, orb]] - \[Mu]) * (LatticeEnergies[[k, orb, orb]] + \[CapitalSigma]0)
		, {k, LE}];
		(* if there is orbital symmetry, avoid the calculation for all the other orbitals *)
		If[OrbitalSymmetry,
			Ekin = Norb * Ekin;
			Break[];
		];
	, {orb, 1, Norb}];
	Chop[Total[Ekin], 10^(-8)]
];

(* EdMode = "Raman" means that the orbitals are not coupled by tunneling terms *)
KineticEnergyRaman[\[Mu]_, LatticeEnergies_, LatticeWeights_, \[CapitalSigma]_, i\[Omega]_, OrbitalSymmetry_, OptionsPattern[]] := Module[
	{flavdist, Energies, Norb, f, LE},
	Norb = Length[LatticeEnergies[[1]]];
	f = If[Norb == 1, Length[\[CapitalSigma][[1]]], Length[\[CapitalSigma][[1,1]]]];
	LE = Length[LatticeWeights];
	(* if flavor distribution has been computed already, use it instead of computing it twice *)
	If[OptionValue["FlavorDistribution"] != {}, flavdist = OptionValue["FlavorDistribution"]];
	Table[
		(* <cdg_k\[Alpha] c_k\[Beta]> for the effective flavors for a given orbital *)
		If[OptionValue["FlavorDistribution"] == {}, 
			flavdist = MomentumDistributedDensityRaman[Range[LE], LatticeEnergies[[All, orb, orb]], \[Mu], If[Norb==1||OrbitalSymmetry, \[CapitalSigma], \[CapitalSigma][[orb]]], i\[Omega]];
		];
		Sum[
			LatticeWeights[[i]] * Tr[ LatticeEnergies[[i, orb, orb]] . Re[ flavdist[[If[OrbitalSymmetry,1,orb], i]] ] ]
		, {i, LE}]
	, {orb, Norb}]
];
Options[KineticEnergyRaman] = {"FlavorDistribution" -> {}};

KineticEnergyInterorbNormal[\[Mu]_, LatticeEnergies_, LatticeWeights_, \[CapitalSigma]_, i\[Omega]_] := Module[
	{Norb, \[CapitalSigma]0, TMats, NMats, LE, G, Ekin},
	Norb = Length[\[CapitalSigma]];
	\[CapitalSigma]0 = Last /@ \[CapitalSigma];
	TMats = Im[ i\[Omega][[2]] - i\[Omega][[1]] ] / (2*Pi);
	NMats = Length[i\[Omega]];
	LE = Length[LatticeWeights];
	(* one factor 2 accounts for spin degeneracy, the other one for sum over negative Matsubara frequencies *)
	Ekin = 2. * 2. * TMats * Sum[
		LatticeWeights[[k]] * Tr[
			(LatticeEnergies[[k]] - \[Mu]) * Sum[
				Re[ Inverse[ (i\[Omega][[n]] + \[Mu])*IdentityMatrix[Norb] - LatticeEnergies[[k]] - DiagonalMatrix[\[CapitalSigma][[All, n]]] ] ]
				- (LatticeEnergies[[k]] + \[CapitalSigma]0)/(i\[Omega][[n]]^2)
		, {n, NMats}]]
	, {k, LE}] + 
	0.5 * Sum[
		LatticeWeights[[k]] * Tr[(LatticeEnergies[[k]] - \[Mu])]
	, {k, LE}] - 
	(1./(2.*TMats)) * Sum[
		LatticeWeights[[k]] * Tr[(LatticeEnergies[[k]] - \[Mu]) * (LatticeEnergies[[k]] + \[CapitalSigma]0)]
	, {k, LE}];
	Chop[Ekin, 10^(-8)]
];

KineticEnergySuperc[\[Mu]_, LatticeEnergies_, LatticeWeights_, \[CapitalSigma]_, i\[Omega]_, OrbitalSymmetry_] := Module[
	{Norb, \[CapitalSigma]0, TMats, NMats, LE, \[Sigma]0 = IdentityMatrix[2], \[Sigma]3 = PauliMatrix[3], Ekin},
	Norb = Length[\[CapitalSigma]];
	Ekin = ConstantArray[0., Norb]; (* initialize Norb orbital contributions to the kinetic energy *)
	TMats = Im[ i\[Omega][[2]] - i\[Omega][[1]] ] / (2*Pi);
	NMats = Length[i\[Omega]];
	LE = Length[LatticeWeights];
	Do[
		\[CapitalSigma]0 = Last[\[CapitalSigma][[orb]]];
		\[CapitalSigma]0[[2,2]] = - Conjugate[ \[CapitalSigma]0[[2,2]] ];
		(* here we do not assume spin degeneracy: we sum spin up and spin down contributions *)
		(* one factor 2 accounts for sum over negative Matsubara frequencies *)
		Ekin[[orb]] = 2. * TMats * Sum[
			LatticeWeights[[k]] * (LatticeEnergies[[k, orb, orb]] - \[Mu]) * Sum[
				Re[ Tr[ TwoByTwoInverse[ i\[Omega][[n]]*\[Sigma]0 + (\[Mu] - LatticeEnergies[[k, orb, orb]])*\[Sigma]3 - \[CapitalSigma][[orb, n]] ] . \[Sigma]3 ] ]
				- (2.*LatticeEnergies[[k, orb, orb]] + Tr[\[CapitalSigma]0])/(i\[Omega][[n]]^2)
			, {n, NMats}]
		, {k, LE}] + 
		0.5 * Sum[
			LatticeWeights[[k]] * (LatticeEnergies[[k, orb, orb]] - \[Mu])
		, {k, LE}] - 
		(1./(4.*TMats)) * Sum[
			LatticeWeights[[k]] * (LatticeEnergies[[k, orb, orb]] - \[Mu]) * (2.*LatticeEnergies[[k, orb, orb]] + Tr[\[CapitalSigma]0])
		, {k, LE}];
		(* if there is orbital symmetry, avoid the calculation for all the other orbitals *)
		If[OrbitalSymmetry,
			Ekin = Norb * Ekin;
			Break[];
		];
	, {orb, 1, Norb}];
	Chop[Total[Ekin], 10^(-8)]
];

KineticEnergyFullSuperc[\[Mu]_, LatticeEnergies_, LatticeWeights_, \[CapitalSigma]_, i\[Omega]_] := Module[
	{Norb, \[CapitalSigma]0, TMats, NMats, Energies, LE, \[Sigma]0 = IdentityMatrix[4], \[Sigma]3 = DiagonalMatrix[{1,-1,1,-1}], Ekin},
	\[CapitalSigma]0 = Last[\[CapitalSigma]];
	Norb = Length[\[CapitalSigma]0]/2;
	TMats = Im[ i\[Omega][[2]] - i\[Omega][[1]] ] / (2*Pi);
	NMats = Length[i\[Omega]];
	LE = Length[LatticeWeights];
	Energies = ConstantArray[0., {LE, 2*Norb, 2*Norb}];
	Do[
		Do[
			Energies[[All, 2*(orb1-1)+1, 2*(orb2-1)+1]] = LatticeEnergies[[All, orb1, orb2]];
			Energies[[All, 2*(orb1-1)+2, 2*(orb2-1)+2]] = - LatticeEnergies[[All, orb1, orb2]];
		, {orb2, Norb}];
	, {orb1, Norb}];
	(**)
	Ekin = 2. * TMats * Sum[
		LatticeWeights[[k]] * Tr[ (Energies[[k]] - \[Mu]*\[Sigma]0) . Sum[
			Re[ Inverse[ i\[Omega][[n]]*\[Sigma]0 + \[Mu]*\[Sigma]3 - Energies[[k]] - \[CapitalSigma][[n]] ] ] 
			- Tr[ (Energies[[k]] + \[CapitalSigma]0) . \[Sigma]3 ]/(i\[Omega][[n]]^2)
		, {n, NMats}] ]
	, {k, LE}] + 
	0.5 * Sum[
		LatticeWeights[[k]] * Tr[ (Energies[[k]] - \[Mu] . \[Sigma]3) . \[Sigma]3 ]
	, {k, LE}] - 
	(1./(4.*TMats)) * Sum[
		LatticeWeights[[k]] * Tr[ (Energies[[k]] - \[Mu]*\[Sigma]3) . (Energies[[k]] + \[CapitalSigma]0) ]
	, {k, LE}];
	Ekin
];

(* call the proper function depending on EdMode *)
KineticEnergy[\[Mu]_, LatticeEnergies_, LatticeWeights_, \[CapitalSigma]_, i\[Omega]_, EdMode_, OptionsPattern[]] := Which[
	EdMode == "Normal", KineticEnergyNormal[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega], OptionValue["OrbitalSymmetry"]],
	EdMode == "Superc", KineticEnergySuperc[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega], OptionValue["OrbitalSymmetry"]],
	EdMode == "Raman", KineticEnergyRaman[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega], OptionValue["OrbitalSymmetry"], "FlavorDistribution" -> OptionValue["FlavorDistribution"]],
	EdMode == "InterorbNormal", KineticEnergyInterorbNormal[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega]],
	EdMode == "FullSuperc", KineticEnergyFullSuperc[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega]]
];
Options[KineticEnergy] = {"OrbitalSymmetry" -> False, "FlavorDistribution" -> {}};


(* flavor current for Raman coupled systems *)
FlavorCurrent[t_, \[Gamma]_, \[Sigma]_, a_, flavordistribution_, LatticeType_, LatticeDim_, NumberOfPoints_, SublatticesQ_] := Module[
	{LE, BZ, MBZ, Iflavor, f, m},
	(* get Brillouin zone *)
	LE = Floor[NumberOfPoints^(1./LatticeDim)];
	If[!SublatticesQ, 
		f = Length[flavordistribution[[1]]];
		BZ = BrillouinZone[LE, LatticeDim, LatticeType];,
	(* else if SublatticeQ *)
		f = Length[flavordistribution[[1]]]/2;
		MBZ = MagneticBrillouinZone[LE, LatticeDim, LatticeType];
	];
	(* spin projection quantum number labeled by \[Sigma]: e.g. with f=2, \[Sigma]=1 -> m=-1/2; \[Sigma]=2 -> m=1/2 *)
	m = -(\[Sigma]-(f+1)/2);
	(* compute flavor current *)
	If[!SublatticesQ,
		Which[
			LatticeType == "Hypercubic",
			Iflavor = (2.*t / (LE^LatticeDim)) * (
				Table[Sin[k[[a]] + m*\[Gamma][[a]]], {k, BZ}] . 
				Re[ flavordistribution[[All, \[Sigma], \[Sigma]]] ]
			);,
		(* ----------------------------------------- *)
			LatticeType == "Bethe",
			Return["Error. Bethe lattice is incompatible with a gauge field. "]
		];,
	(* else, if SublatticesQ *)
		Which[
			LatticeType == "Hypercubic" && LatticeDim == 1,
			Iflavor = -2.0*(t/(LE^LatticeDim)) * Re[ Table[
				- I*Exp[I*m*\[Gamma][[1]]] + I*Exp[-I*m*\[Gamma][[1]]]*Exp[-2.0*I*k[[1]]]
			, {k, MBZ}] . flavordistribution[[All, \[Sigma], 2+\[Sigma]]] ],
		(* ----------------------------------------- *)
			LatticeType == "Hypercubic" && LatticeDim == 2,
			Iflavor = -2.0*(t/(LE^LatticeDim)) * Re[ Table[
				If[a == 1, (* x direction *)
					- I * Exp[I*m*\[Gamma][[1]]] 
					+ I * Exp[-I*m*\[Gamma][[1]]]*Exp[-2I*k[[1]]],
				(* else, y direction *)
					- I * Exp[I*m*\[Gamma][[2]]]*Exp[I*(-k[[1]]+k[[2]])] 
					+ I * Exp[-I*m*\[Gamma][[2]]]*Exp[-I*(k[[1]]+k[[2]])]
				]
			, {k, MBZ}] . flavordistribution[[All, \[Sigma], 2+\[Sigma]]] ],
		(* ----------------------------------------- *)
			LatticeType == "Bethe",
			Return["Error. Bethe lattice is incompatible with a gauge field. "]
		];
	];
	Iflavor
];

(* non interacting chiral current in thermodynamic limit, d=1 *)
ChiralCurrentNonInteracting[t_, \[CapitalOmega]_, \[Gamma]_]:= 
If[
	\[CapitalOmega] <= 2.0*t*Cos[\[Gamma]/2], (* metal *)
		2.0*(-0.3183098861837907`)*(-((2.` Cos[\[Gamma]/2]^2 Sin[\[Gamma]/2])/Sqrt[1-0.25` \[CapitalOmega]^2-Sin[\[Gamma]/2]^2])-2.` Sin[\[Gamma]/2] Sqrt[1-0.25` \[CapitalOmega]^2-Sin[\[Gamma]/2]^2]+\[CapitalOmega] (0.5` Cot[\[Gamma]/2] (EllipticE[\[Pi]-ArcCos[Sqrt[0.25` \[CapitalOmega]^2+Sin[\[Gamma]/2]^2]],-((4.` Sin[\[Gamma]/2]^2)/\[CapitalOmega]^2)]-EllipticF[\[Pi]-ArcCos[Sqrt[0.25` \[CapitalOmega]^2+Sin[\[Gamma]/2]^2]],-((4.` Sin[\[Gamma]/2]^2)/\[CapitalOmega]^2)])-0.5` Cot[\[Gamma]/2] (EllipticE[ArcCos[Sqrt[0.25` \[CapitalOmega]^2+Sin[\[Gamma]/2]^2]],-((4.` Sin[\[Gamma]/2]^2)/\[CapitalOmega]^2)]-EllipticF[ArcCos[Sqrt[0.25` \[CapitalOmega]^2+Sin[\[Gamma]/2]^2]],-((4.` Sin[\[Gamma]/2]^2)/\[CapitalOmega]^2)])+(Cos[\[Gamma]/2] Sin[\[Gamma]/2] Sqrt[1+(4.` Sin[\[Gamma]/2]^2 (1-0.25` \[CapitalOmega]^2-Sin[\[Gamma]/2]^2))/\[CapitalOmega]^2])/(Sqrt[1-0.25` \[CapitalOmega]^2-Sin[\[Gamma]/2]^2] Sqrt[0.25` \[CapitalOmega]^2+Sin[\[Gamma]/2]^2]))),
	(* else, insulator *)
		2.0*(-0.15915494309189535`) \[CapitalOmega] Cot[\[Gamma]/2] (EllipticE[\[Pi],-((4.` Sin[\[Gamma]/2]^2)/\[CapitalOmega]^2)]-EllipticF[\[Pi],-((4.` Sin[\[Gamma]/2]^2)/\[CapitalOmega]^2)])
];


End[]

EndPackage[]
