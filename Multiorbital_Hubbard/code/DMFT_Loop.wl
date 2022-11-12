(* ::Package:: *)

(*                   DMFT LOOP                     *)
Do[
	ClearAll[Hsectors, EgsSectorList, GsSectorList];
	
	(* First print *)
	Print[Style["DMFT Loop n. ", 20, Bold, Red], Style[DMFTiterator, 20, Bold, Red]];
	Print[Style["\t\t Exact Diagonalization start", 16, Bold, Orange]];
	Print["Bath parameters: ", BathParameters];

	(* Build and diagonalize the AIM Hamiltonian + print timing *)
	Print["E.D. time: ", AbsoluteTiming[
	
		Hsectors = HImp[Norb, HnonlocBlocks, HlocBlocks, BathParameters, InteractionParameters, EdMode];
		eigs = Map[Eigs[#, "Temperature" -> T, "MinLanczosDim"->100, "MaxIterations" -> 2000]&, Hsectors];
		EgsSectorList = eigs[[All, 1]];
		GsSectorList = eigs[[All, 2]];
		ClearAll[eigs];
		
	]," sec.\n"];
	
	
	(* -------------------------- *)
	(* ZERO TEMPERATURE CALCULATIONS *)
	(* -------------------------- *)
	If[T == 0,
		
		(* COMPUTE THE GROUND STATE *)
		Print["Computing ground state and Green functions..."];
		(* ground state energy (lowest of all the sectors) *)
		Egs = Min[Flatten[EgsSectorList]];
		(* sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy *)
		GsSectorIndex = Position[
			EgsSectorList,
			_?((Abs[# - Egs] < DegeneracyThreshold)&)
		];
		(* list of all the degenerate ground states *)
		Gs = MapApply[GsSectorList[[##]]&, GsSectorIndex];
		(* list of quantum numbers of the degenerate ground states *)
		GsQns = QnsSectorList[[GsSectorIndex[[All, 1]]]];	
		
		(* print relevant information about the ground state *)
		Print["\t\t Ground state info:\n", "Egs = ", Egs, "   Quantum numbers = ", GsQns];
		(* print some observables *)
		Do[
			Print["impurity density [orb="<>ToString[orb]<>"] = ", Sum[Density[L, f, Norb, 1, \[Sigma], orb, Sectors, EgsSectorList, GsSectorList, T], {\[Sigma], f}] ];
			Print["impurity double occupancy [orb="<>ToString[orb]<>"] = ", SquareDensity[L, f, Norb, {1,1}, {1,2}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T] ];
			If[EdMode == "Superc", Print["order parameter [orb="<>ToString[orb]<>"] = ", CdgCdg[L, f, Norb, {1,1}, {1,2}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T] ];]
		, {orb, Norb}];
		
		
		(* --------------------------------------- *)
		(* if there is orbital symmetry, you compute many body functions just for ONE representative orbital *)
		(* --------------------------------------- *)
		Which[
			(EdMode == "Normal" || EdMode == "Superc") && OrbitalSymmetry,
			(* identify independent parameters, i.e. the minimal set of bath parameters that you need to compute stuff *)
			IndependentParameters = TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode];
			(* G^-1(i\[Omega]) *)
			InverseG = Mean[MapApply[
				InverseGreenFunction[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
				{Gs, GsQns}\[Transpose]
			]];
			(* G_0^-1(i\[Omega]) *)
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0]; (* memorize the previous one *)
			InverseG0 = (Weiss/.Thread[symbols -> IndependentParameters])/.{z -> #}&/@i\[Omega];
			(* \[CapitalSigma](i\[Omega]) *)
			\[CapitalSigma] = InverseG0 - InverseG;
			(* G_loc(i\[Omega]) *)
			LocalG = LocalGreenFunction[LatticeEnergies[[All,1,1]], LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode];
			
			(* print stuff *)
			Print["Quasiparticle weight z = ", QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode] ];
			Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
			Print["S.C. time: ", First@AbsoluteTiming[
			
			(* Self consistency *)
			NewBathParameters = ReshapeBathParameters[L, f, Norb,	
				SelfConsistency[W[[1]], \[Mu], Weiss, symbols, z, IndependentParameters, LocalG, \[CapitalSigma], i\[Omega], EdMode, 
				Lattice -> LatticeType, LatticeDimension -> LatticeDim, Minimum -> "Local", NumberOfFrequencies -> 500, MaxIterations -> 2000, AccuracyGoal -> 7],
			OrbitalSymmetry, EdMode];
			
			], " sec." ];
			
			(* compute error *)
			error = DMFTError[InverseG0[[;;2000]], InverseG0old[[;;2000]], EdMode];,
		
		(* -------------------- *)	
		(* if NO ORBITAL SYMMETRY *)
		(* -------------------- *)
			(EdMode == "Normal" || EdMode == "Superc") && !OrbitalSymmetry,
			IndependentParameters = Table[
				TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode],
			{orb, Norb}];
			(* { G^-1(i\[Omega])_orb=1 , G^-1(i\[Omega])_orb=2 , ...} *)
			InverseG = Table[
				Mean[MapApply[
					InverseGreenFunction[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
					{Gs, GsQns}\[Transpose]
				]], {orb, Norb}];
			(* { Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=1] , Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=2] , ...} *)
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0];
			InverseG0 = Table[
				(Weiss/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode]])/.{z -> #}&/@i\[Omega],
				{orb, Norb}];
			(* { \[CapitalSigma]Subscript[(i\[Omega]), orb=1] , \[CapitalSigma]Subscript[(i\[Omega]), orb=2] , ...} *)
			\[CapitalSigma] = InverseG0 - InverseG;
			(* { Subscript[G, loc]Subscript[(i\[Omega]), orb=1] , Subscript[G, loc]Subscript[(i\[Omega]), orb=2] , ...} *)
			LocalG = Table[
				LocalGreenFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma][[orb]], i\[Omega], EdMode]
			, {orb, Norb}];
			(* Self consistency *)
			Print["Quasiparticle weight z = ", Table[QuasiparticleWeight[\[CapitalSigma][[orb]], i\[Omega], EdMode], {orb, Norb}] ];
			Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
			Print["S.C. time: ", First@AbsoluteTiming[
			
			NewBathParameters = ReshapeBathParameters[L, f, Norb, Table[
				SelfConsistency[W[[orb]], \[Mu], Weiss, symbols, z, IndependentParameters[[orb]], LocalG[[orb]], \[CapitalSigma][[orb]], i\[Omega], EdMode,
				Lattice -> LatticeType, LatticeDimension -> LatticeDim, Minimum -> "Local", NumberOfFrequencies -> 500, MaxIterations -> 2000, AccuracyGoal -> 7]
			, {orb, Norb}], OrbitalSymmetry, EdMode];
			
			], " sec." ];
			
			error = (1./Norb)*Sum[
				DMFTError[InverseG0[[orb]], InverseG0old[[orb]], EdMode],
			{orb, Norb}];
		];
	];
	
	(* -----------------------------*)
	(* FINITE TEMPERATURE CALCULATIONS *)
	(* -----------------------------*)
	If[T != 0,
		Print["not supported."];
		Break[];
	];
	
	(* update bath parameters *)
	BathParameters = Mixing * BathParameters + (1 - Mixing) * NewBathParameters;
	Print["DMFT error: ", ScientificForm[error]];
	AppendTo[ErrorList, error];
	
	(*Print*)
	Print[Style["\t\t Self Consistency completed", 16, Bold, Magenta]];
	Print["----------------------------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------------------------"];

	(* Exit DMFT Loop if convengerce is reached *)
	If[DMFTiterator > DMFTMinIterations && error < DMFTerror && LastIteration, Break[];];
	If[error < DMFTerror, LastIteration = True, (*else*) LastIteration = False];


, {DMFTiterator, DMFTMaxIterations}]
