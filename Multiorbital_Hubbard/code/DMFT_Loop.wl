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
		Gs = MapApply[GsSectorList[[##]]&, GsSectorIndex];
		(* list of quantum numbers of the degenerate ground states *)
		GsQns = QnsSectorList[[GsSectorIndex[[All, 1]]]];	
		
		(* print relevant information about the ground state *)
		Print["\t\t Ground state info:\n", "Egs = ", Egs, "   Quantum numbers = ", GsQns];
		(* print some observables *)
		Print["\n\t\t Observables:"];
		Do[
			Print["Impurity density [orb="<>ToString[orb]<>"] = ", Sum[Density[L, f, Norb, 1, \[Sigma], orb, Sectors, EgsSectorList, GsSectorList, T], {\[Sigma], f}] ];
			Print["Impurity double occupancy [orb="<>ToString[orb]<>"] = ", SquareDensity[L, f, Norb, {1,1}, {1,2}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T] ];
			If[EdMode == "Superc" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc", 
				Print["Order parameter [orb="<>ToString[orb]<>"] = ", CdgCdg[L, f, Norb, {1,1}, {1,2}, {orb,orb}, Sectors, EgsSectorList, GsSectorList, T] ];
			];
		, {orb, Norb}];
		
		If[EdMode == "InterorbSuperc" || EdMode == "FullSuperc", 
				Print["Order parameter interorbital = ", 
				0.5 * (CdgCdg[L, f, Norb, {1,1}, {1,2}, {1,2}, Sectors, EgsSectorList, GsSectorList, T]
				- CdgCdg[L, f, Norb, {1,1}, {1,2}, {2,1}, Sectors, EgsSectorList, GsSectorList, T])];
				(* the minus sign is WRONG. it should come from the second CdgCdg, not put explicitly here. *)
			];
		
		
		Which[
		(* --------------------------------------- *)
		(* if there is orbital symmetry, you compute many body functions just for ONE representative orbital *)
		(* --------------------------------------- *)
			(EdMode == "Normal" || EdMode == "Superc") && OrbitalSymmetry,
			Print["\n Green functions calculation time: ", First@AbsoluteTiming[
			
			(* identify independent parameters, i.e. the minimal set of bath parameters that you need to compute stuff *)
			IndependentParameters = TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode];
			(* G(i\[Omega]) *)
			Gimp = Mean[MapApply[
				GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
				{Gs, GsQns}\[Transpose]
			]];
			(* G^-1(i\[Omega]) *)
			InverseG = InverseGreenFunction[Gimp, EdMode];
			(* G_0^-1(i\[Omega]) *)
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0]; (* memorize the previous one *)
			InverseG0 = ((Weiss/.{\[Mu]eff -> \[Mu] - \[Delta][[1]]})/.Thread[symbols -> IndependentParameters])/.{z -> #}&/@i\[Omega];
			(* \[CapitalSigma](i\[Omega]) *)
			\[CapitalSigma] = InverseG0 - InverseG;
			(* G_loc(i\[Omega]) *)
			LocalGold = If[DMFTiterator == 1, 0*InverseG, LocalG];
			LocalG = LocalGreenFunction[LatticeEnergies[[All,1,1]], LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode];
		
		(* IN PROGRESS ... *)
			If[DMFTiterator > 2,
				LocalG = (1.0 - Mixing) * LocalG + Mixing * LocalGold;
			];
		(* END IN PROGRESS *)	
		
			(* compute and print z *)
			Z = QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode, FitCutoff -> 50];
			Print["Quasiparticle weight z = ", Z];
			
			], " sec."];
			
			
			(* Self consistency *)
			Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
			Print["S.C. time: ", First@AbsoluteTiming[
			
			NewBathParameters = ReshapeBathParameters[L, f, Norb,	
				SelfConsistency[
					W[[1]], \[Mu] - \[Delta][[1]], Weiss/.{\[Mu]eff -> \[Mu] - \[Delta][[1]]}, symbols, z, IndependentParameters, LocalG, \[CapitalSigma], i\[Omega], EdMode, 
					Lattice -> LatticeType, 
					LatticeDimension -> LatticeDim, 
					Minimum -> MinimizationType, 
					Method -> MinimizationMethod,
					NumberOfFrequencies -> CGNMatsubara, 
					MaxIterations -> CGMaxIterations, 
					AccuracyGoal -> CGAccuracy,
					FitWeight -> CGWeight],
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
				Mean[MapApply[
					GreenFunctionImpurity[L, f, Norb, 1, orb, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
					{Gs, GsQns}\[Transpose]
				]], {orb, Norb}];
			(* G^-1(i\[Omega])_orb=1 , G^-1(i\[Omega])_orb=2 ... *)
			InverseG = InverseGreenFunction[#, EdMode] &/@ Gimp;
			(* { Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=1] , Subscript[G, 0]^-1Subscript[(i\[Omega]), orb=2] , ...} *)
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0];
			InverseG0 = Table[(
				(Weiss/.{\[Mu]eff -> \[Mu] - \[Delta][[orb]]})/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, orb, BathParameters, EdMode]])/.{z -> #}&/@i\[Omega]
			, {orb, Norb}];
			(* { \[CapitalSigma]Subscript[(i\[Omega]), orb=1] , \[CapitalSigma]Subscript[(i\[Omega]), orb=2] , ...} *)
			\[CapitalSigma] = InverseG0 - InverseG;
			(* { Subscript[G, loc]Subscript[(i\[Omega]), orb=1] , Subscript[G, loc]Subscript[(i\[Omega]), orb=2] , ...} *)
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
			
			NewBathParameters = ReshapeBathParameters[L, f, Norb, Table[
				SelfConsistency[
					W[[orb]], \[Mu] - \[Delta][[orb]], Weiss/.{\[Mu]eff -> \[Mu] - \[Delta][[orb]]}, symbols, z, IndependentParameters[[orb]], LocalG[[orb]], \[CapitalSigma][[orb]], i\[Omega], EdMode,
					Lattice -> LatticeType, 
					LatticeDimension -> LatticeDim, 
					Minimum -> MinimizationType, 
					Method -> MinimizationMethod,
					NumberOfFrequencies -> CGNMatsubara, 
					MaxIterations -> CGMaxIterations, 
					AccuracyGoal -> CGAccuracy,
					FitWeight -> CGWeight]
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
			Gimp = Mean[MapApply[
				GreenFunctionImpurity[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
				{Gs, GsQns}\[Transpose]
			]];
			InverseG = InverseGreenFunction[Gimp, EdMode];
			InverseG0old = If[DMFTiterator == 1, 0*InverseG, InverseG0];
			InverseG0 = (
				(Weiss/.{\[Mu]eff -> \[Mu] - \[Delta]})/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode]]
			)/.{z -> #} &/@ i\[Omega];
			\[CapitalSigma] = InverseG0 - InverseG;
			LocalG = LocalGreenFunction[LatticeEnergies, LatticeWeights, \[Mu], \[CapitalSigma], i\[Omega], EdMode];
			
			Z = Table[QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode, Orb -> orb], {orb, Norb}];
			Print["Quasiparticle weight z = ", Z];
			
			], " sec."];
			
			(* Self consistency *)
			(* Print["Quasiparticle weight z = ", QuasiparticleWeight[\[CapitalSigma], i\[Omega], EdMode] ]; *)
			Print[Style["\t\t Self Consistency start", 16, Bold, Magenta]];
			Print["S.C. time: ", First@AbsoluteTiming[
			
			NewBathParameters = ReshapeBathParameters[L, f, Norb, 
				SelfConsistency[
					W, \[Mu] - \[Delta], Weiss/.{\[Mu]eff -> \[Mu] - \[Delta]}, symbols, z, IndependentParameters, LocalG, \[CapitalSigma], i\[Omega], EdMode,
					Lattice -> LatticeType, 
					LatticeDimension -> LatticeDim, 
					Minimum -> MinimizationType, 
					Method -> MinimizationMethod,
					NumberOfFrequencies -> CGNMatsubara, 
					MaxIterations -> CGMaxIterations, 
					AccuracyGoal -> CGAccuracy,
					FitWeight -> CGWeight]
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
	
	(* update bath parameters *)
	BathParameters = NewBathParameters;
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
