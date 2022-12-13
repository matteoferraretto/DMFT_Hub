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
		eigs = Map[
			Eigs[#, 
				"Temperature" -> T, 
				"MinLanczosDim" -> MinLanczosDim, 
				"DegeneracyThreshold" -> DegeneracyThreshold,
				"MaxIterations" -> MaxLanczosIter, 
				"MinEigenvalues" -> MinNumberOfEigs
			]&, Hsectors];
		EgsSectorList = eigs[[All, 1]];
		GsSectorList = eigs[[All, 2]];
		ClearAll[eigs];
		
	]," sec.\n"];
	
	
	(* -------------------------- *)
	(* ZERO TEMPERATURE CALCULATIONS *)
	(* -------------------------- *)
	If[T == 0,
		
		(* COMPUTE THE GROUND STATE *)
		(* ground state energy (lowest of all the sectors) *)
		Egs = Min[Flatten[EgsSectorList]];
		(* sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy *)
		GsSectorIndex = Position[
			EgsSectorList,
			_?((Abs[# - Egs] < DegeneracyThreshold)&)
		];
		(* list of all the degenerate ground states *)
		Gs = Apply[GsSectorList[[##]]&, GsSectorIndex, {1}];
		(* list of quantum numbers of the degenerate ground states *)
		GsQns = QnsSectorList[[GsSectorIndex[[All, 1]]]];	
		
		(* print relevant information about the ground state *)
		Print["\t\t Ground state info:\n", "Egs = ", Egs, "   Quantum numbers = ", GsQns];
		(* print some observables *)
		Print["\n\t\t Observables:"];
		density = Table[ Sum[Density[L, f, Norb, 1, \[Sigma], orb, Sectors, EgsSectorList, GsSectorList, T], {\[Sigma], f}], {orb, Norb}];
		Print["Impurity density = ", density ];
		docc = Table[ SquareDensity[L, f, Norb, {1,1}, {1,2}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T], {orb, Norb}];
		Print["Impurity double occupancy = ", docc ];
		If[EdMode == "Superc" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
			\[Phi] = Table[ CdgCdg[L, f, Norb, {1,1}, {1,2}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T], {orb, Norb}];
			Print["Order parameter = ", \[Phi] ];
			If[EdMode != "Superc",
				\[CapitalXi] = {
					CdgCdg[L, f, Norb, {1,1}, {1,2}, {1,2}, Sectors, EgsSectorList, GsSectorList, T],
					CdgCdg[L, f, Norb, {1,1}, {1,2}, {2,1}, Sectors, EgsSectorList, GsSectorList, T]
				};
				Print["Order parameter interorbital {<adg_up bdg_dw>, <bdg_up adg_dw>} = ", \[CapitalXi]];
			]
		];
		
		
		Which[
		(* ----------------------------------------------------------------------------------------- *)
		(* if there is orbital symmetry, you compute many body functions just for ONE representative orbital *)
		(* ----------------------------------------------------------------------------------------- *)
			(EdMode == "Normal" || EdMode == "Superc") && OrbitalSymmetry,
			Print["\n Green functions calculation time: ", First@AbsoluteTiming[
			
			(* identify independent parameters, i.e. the minimal set of bath parameters that you need to compute stuff *)
			IndependentParameters = TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode];
			(* G(i\[Omega]) *)
			Gimp = Mean[Apply[
				GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
				{Gs, GsQns}\[Transpose]
			, {1}]];
			(* G^-1(i\[Omega]) *)
			InverseG = InverseGreenFunction[Gimp, EdMode];
			(* G_0^-1(i\[Omega]) *)
			If[DMFTiterator == 1, Weiss = WeissField[L, f, Norb, \[Mu] - \[Delta][[1]], symbols, z, EdMode]; ];
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0]; (* memorize the previous one *)
			InverseG0 = (Weiss/.Thread[symbols -> IndependentParameters])/.{z -> #} &/@ i\[Omega];
			(* \[CapitalSigma](i\[Omega]) *)
			\[CapitalSigma]old = If[DMFTiterator == 1, 0*InverseG, \[CapitalSigma]];
			\[CapitalSigma] = InverseG0 - InverseG;
			(* G_loc(i\[Omega]) *)
			LocalGold = If[DMFTiterator == 1, 0*InverseG, LocalG];
			LocalG = LocalGreenFunction[LatticeEnergies[[All,1,1]], LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode];
		
			(* compute and print z *)
			Z = QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode, FitCutoff -> 50];
			Print["Quasiparticle weight z = ", Z];
			
			], " sec."];
			
			
			(* Self consistency *)
			Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
			Print["S.C. time: ", First@AbsoluteTiming[
			
			BathParameters = ReshapeBathParameters[L, f, Norb,	
				SelfConsistency[
					W[[1]], \[Mu] - \[Delta][[1]], Weiss, symbols, z, IndependentParameters, 
					LocalG, LocalGold, \[CapitalSigma], \[CapitalSigma]old, i\[Omega], EdMode,
					Mix -> If[DMFTiterator > 2, Mixing, (* else *) 0.0],
					Lattice -> LatticeType, 
					LatticeDimension -> LatticeDim, 
					Minimum -> MinimizationType, 
					Method -> MinimizationMethod,
					NumberOfFrequencies -> CGNMatsubara, 
					MaxIterations -> CGMaxIterations, 
					AccuracyGoal -> CGAccuracy,
					FitWeight -> CGWeight
				],
			OrbitalSymmetry, EdMode];
			
			], " sec." ];
			
			(* compute error *)
			error = DMFTError[InverseG0[[;;CGNMatsubara]], InverseG0old[[;;CGNMatsubara]], EdMode];,
		
		(* -------------------- *)	
		(* if NO ORBITAL SYMMETRY *)
		(* -------------------- *)
			(EdMode == "Normal" || EdMode == "Superc") && !OrbitalSymmetry,
			Print["\n Green functions calculation time: ", First@AbsoluteTiming[
			
			IndependentParameters = Table[
				TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode],
			{orb, Norb}];
			(* { G(i\[Omega])_orb=1 , G(i\[Omega])_orb=2 , ...} *)
			Gimp = Table[
				Mean[Apply[
					GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
					{Gs, GsQns}\[Transpose]
				, {1}]], {orb, Norb}];
			(* G^-1(i\[Omega])_orb=1 , G^-1(i\[Omega])_orb=2 ... *)
			InverseG = InverseGreenFunction[#, EdMode] &/@ Gimp;
			(* { Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=1] , Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=2] , ...} *)
			If[DMFTiterator == 1, 
				Weiss = Table[ WeissField[L, f, Norb, \[Mu] - \[Delta][[orb]], symbols, z, EdMode], {orb, Norb}]
			];
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0];
			InverseG0 = Table[(
				Weiss[[orb]]/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode]])/.{z -> #}&/@i\[Omega]
			, {orb, Norb}];
			(* { \[CapitalSigma]Subscript[(i\[Omega]), orb=1] , \[CapitalSigma]Subscript[(i\[Omega]), orb=2] , ...} *)
			\[CapitalSigma]old = If[DMFTiterator == 1, 0*InverseG, \[CapitalSigma]];
			\[CapitalSigma] = InverseG0 - InverseG;
			(* { Subscript[G, loc]Subscript[(i\[Omega]), orb=1] , Subscript[G, loc]Subscript[(i\[Omega]), orb=2] , ...} *)
			LocalGold = If[DMFTiterator == 1, 0*InverseG, LocalG];
			LocalG = Table[
				LocalGreenFunction[LatticeEnergies[[All, orb, orb]], LatticeWeights, \[Mu], \[CapitalSigma][[orb]], i\[Omega], EdMode]
			, {orb, Norb}];
			
			(* z *)
			Z = Table[ QuasiparticleWeight[\[CapitalSigma][[orb]], i\[Omega], EdMode, FitCutoff -> 50], {orb, Norb} ];
			Print["Quasiparticle weight z = ", Z];
			
			], " sec."];
			
			
			(* Self consistency *)
			Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
			Print["S.C. time: ", First@AbsoluteTiming[
			
			BathParameters = ReshapeBathParameters[L, f, Norb, Table[
				SelfConsistency[
					W[[orb]], \[Mu] - \[Delta][[orb]], Weiss[[orb]], symbols, z, IndependentParameters[[orb]], 
					LocalG[[orb]], LocalGold[[orb]], \[CapitalSigma][[orb]], \[CapitalSigma]old[[orb]], i\[Omega], EdMode,
					Mix -> If[DMFTiterator > 2, Mixing, (* else *) 0.0],
					Lattice -> LatticeType, 
					LatticeDimension -> LatticeDim, 
					Minimum -> MinimizationType, 
					Method -> MinimizationMethod,
					NumberOfFrequencies -> CGNMatsubara, 
					MaxIterations -> CGMaxIterations, 
					AccuracyGoal -> CGAccuracy,
					FitWeight -> CGWeight
				]
			, {orb, Norb}], OrbitalSymmetry, EdMode];
			
			], " sec." ];
			
			error = (1./Norb)*Sum[
				DMFTError[InverseG0[[orb]], InverseG0old[[orb]], EdMode],
			{orb, Norb}];,
		
		(* ------------------ *)
		(*   Interorb Superc   *)
		(* ------------------ *)
			(EdMode == "InterorbSuperc" || EdMode == "FullSuperc") && !OrbitalSymmetry,
			Print["\n Green functions calculation time: ", First@AbsoluteTiming[
			
			IndependentParameters = TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode];
			(* G(i\[Omega]) *)
			Gimp = Mean[Apply[
				GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
				{Gs, GsQns}\[Transpose]
			, {1}]];
			InverseG = InverseGreenFunction[Gimp, EdMode];
			(* G0^-1(i\[Omega]) *)
			If[DMFTiterator == 1,
				Weiss = WeissField[L, f, Norb, \[Mu] - \[Delta], symbols, z, EdMode]
			];
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0];
			InverseG0 = (
				Weiss/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode]]
			)/.{z -> #} &/@ i\[Omega];
			(* \[CapitalSigma](i\[Omega]) *)
			\[CapitalSigma]old = If[DMFTiterator == 1, 0*InverseG, \[CapitalSigma]];
			\[CapitalSigma] = InverseG0 - InverseG;
			(* Gloc(i\[Omega]) *)
			LocalGold = If[DMFTiterator == 1, 0*InverseG, LocalG];
			LocalG = LocalGreenFunction[LatticeEnergies, LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode];
			
			Z = Table[QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode, Orb -> orb], {orb, Norb}];
			Print["Quasiparticle weight z = ", Z];
			
			], " sec."];
			
			(* Self consistency *)
			(* Print["Quasiparticle weight z = ", QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode] ]; *)
			Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
			Print["S.C. time: ", First@AbsoluteTiming[
			
			BathParameters = ReshapeBathParameters[L, f, Norb, 
				SelfConsistency[
					W[[1]], 
					\[Mu] - \[Delta], 
					Weiss, 
					symbols, z, IndependentParameters, 
					LocalG, LocalGold, \[CapitalSigma], \[CapitalSigma]old, i\[Omega], EdMode,
					Mix -> If[DMFTiterator > 2, Mixing, (* else *) 0.0],
					Lattice -> LatticeType, 
					LatticeDimension -> LatticeDim, 
					Minimum -> MinimizationType, 
					Method -> MinimizationMethod,
					NumberOfFrequencies -> CGNMatsubara, 
					MaxIterations -> CGMaxIterations, 
					AccuracyGoal -> CGAccuracy,
					FitWeight -> CGWeight
				]
			, OrbitalSymmetry, EdMode];
			
			error = DMFTError[InverseG0, InverseG0old, EdMode];
			
			], " sec." ];
		];
	];
	
	(* -----------------------------*)
	(* FINITE TEMPERATURE CALCULATIONS *)
	(* -----------------------------*)
	If[T != 0,
		Print["not supported."];
		Break[];
	];
	
	(* store new bath parameters and error *)
	WriteOutput[True, OutputDirectory, "hamiltonian_restart", BathParameters];
	AppendTo[ErrorList, error];
	WriteOutput[True, OutputDirectory, "error", ErrorList];
	Print["DMFT error: ", ScientificForm[error]];
	
	(*Print*)
	Print[Style["\t\t Self Consistency completed", 16, Bold, Magenta]];
	Print["----------------------------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------------------------"];

	(* Exit DMFT Loop if convengerce is reached *)
	If[DMFTiterator > DMFTMinIterations && error < DMFTerror && LastIteration, Converged = True; Break[];];
	If[error < DMFTerror, LastIteration = True, (*else*) LastIteration = False];

, {DMFTiterator, DMFTMaxIterations}]
