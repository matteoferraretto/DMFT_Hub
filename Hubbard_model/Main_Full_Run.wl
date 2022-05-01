(* ::Package:: *)

(* ::Subtitle:: *)
(*LOAD THE DMFT PACKAGE*)


(* ::Input::Initialization:: *)
FolderPath=NotebookDirectory[];
<<(FolderPath<>"Package.wl");
?"DMFT`*"


(* ::Subtitle:: *)
(*INPUT PARAMETERS FOR DMFT LOOP*)


(* ::Input::Initialization:: *)
(*             GENERAL INPUT              *)
Nbath=4;(*number of bath sites*)
Nimp=1;(*number of impurity sites*)
L=Nimp+Nbath;(*total number of sites: bath+impurity*)
Nf=2;(*number of flavours*)
EdMode="Superc";(* "Normal" = no symmetry breaking;  "Superc" = bath exchanging pairs with reservoir *)

(*      INPUT PHYSICAL PARAMETERS        *)
Umin = 0.2;	Umax = 5.00;	dU = 0.2;
Ulist = Table[-U,{U,Umin,Umax,dU}];(* default protocol of U values *)
\[CapitalDelta]0=0.5;(* starting value of symmetry breaking field if EdMode="Superc" *)
\[Mu]=0;(*chemical potential*)

(*         INFO ON DMFT LOOPS         *)
DMFTerror=10^(-5);(*convergence threshold*)
DMFTMinIterations=1;(*minimum number of iterations*)
DMFTMaxIterations=100;(*maximum number of iterations*)
DegeneracyThreshold=10^(-9);(*two states are degenerate when |Subscript[E, 1]-Subscript[E, 2]|<DegeneracyThreshold *)
Mixing=0;(*mixing parameter: at iteration n set {Subscript[e, n],Subscript[V, n]}=Mixing*{Subscript[e, n-1],Subscript[V, n-1]}+(1-Mixing)*{Subscript[e, n],Subscript[V, n]} *)
LFit=1000;(*number of Matsubara frequencies used for the \[Chi]^2 fit*)
LE=500;(* number of channels used to divide the interval [-DBethe,+DBethe] for the rough calculation of \[Integral]d\[Epsilon]D[\[Epsilon]]... *)

(*      INFO ON MATSUBARA AND REAL FREQUENCIES         *)
T=0.001;(*fictitious temperature to define Matsubara frequencies*)
NMatsubara=5000;(*Total number of Matsubara frequencies*)
i\[Omega]=Table[(2*n+1)*Pi*I*T,{n,0,NMatsubara-1}];(*list of Matsubara frequencies*)
\[Omega]min=-5.;\[Omega]max=5.;(*set min and max value for the set of real frequencies*)
NReal=10000;(*number of real frequencies*)
d\[Omega]=(\[Omega]max-\[Omega]min)/NReal;(*real frequency step*)
\[Omega]=Table[\[Omega]min+n*d\[Omega],{n,0,NReal}];(*list of real frequencies*)
\[Eta]=0.125;(*small shift of the pole in the imaginary axis: this avoids singularities, but introduces an artificial broadening of the spectrum*)

(*          INFO ON OUTPUT STORAGE FILE PATH         *)
LoadHamiltonianQ = False; (*set True if you want to load the Hamiltonian from a file, False if you want to generate it on the fly and then store it on a file *)
ImpHBlocksFile = FolderPath<>"ImpHBlocks_L="<>ToString[L]<>".mx";
ImpHLocalFile = FolderPath<>"ImpHLocal_L="<>ToString[L]<>".mx";
StoreObservablesQ = True;



(* ::Subtitle:: *)
(*DMFT LOOP*)


(* GET SECTORS *)
QnsSectorList = SectorList[L, EdMode];(* list of quantum numbers of all the sectors {n,nup} or sz *)
DimSectorList = DimSector[L, #, EdMode]&/@QnsSectorList;(*list of dimensions of all the sectors*)

(* GET HAMILTONIANS *)
{impHblocks,impHlocal} = GetHamiltonian[L, Nf, QnsSectorList, LoadHamiltonianQ, ImpHBlocksFile, ImpHLocalFile, EdMode];

(* INPUT RECAP *)
Print[Style["Recap of input:",16,Bold]];
Print["Nbath: ",Nbath,". Nsectors: ",Dim@QnsSectorList,". Dim. of the largest sector: ",Max@DimSectorList];
NU = Dim@Ulist;
Print["Number of U values in the protocol: ",NU,"\n"];

(* LOOP OVER THE PROTOCOL VALUES OF U *)
Do[
	
	U = Ulist[[u]];
	LastIteration=False;(*allows to do one more iteration after convergence threshold is reached*)
	Converged=False;(*True if DMFT has converged, false otherwise*)
	ErrorList={};(*list of DMFT errors*) 

	(* GET BATH PARAMETERS *)
	InitializeBathMode = If[u == 1, "Default", (*else*) FolderPath<>"hamiltonian_restart_U="<>ToString[Abs@Ulist[[u-1]]]<>".txt" ];
	Which[
		EdMode=="Normal",
		{e,V} = StartingBath[InitializeBathMode, Nbath, EdMode];
		Parameters = Flatten@{e,V},
	(* ----------------------------------------------------------------- *)
		EdMode=="Superc", 
		{e,V,\[CapitalDelta]} = StartingBath[InitializeBathMode, Nbath, EdMode]; 
		\[CapitalDelta]=\[CapitalDelta]0*\[CapitalDelta]; 
		Parameters=Flatten@{e,V,\[CapitalDelta]}
	];
	Nparams = Dim@Parameters;(* total number of parameters *)

	(*                   DMFT LOOP                     *)
	Do[

		(* First print *)
		Print[
			Style["DMFT Loop n. ",20,Bold,Red],
			Style[DMFTiterator,20,Bold,Red],
			Style[". U = ",20,Bold,Red],
			Style[U,20,Bold,Red]	
		];
		Print["----------------------------------------------------------------------------------------"];
		Print[Style["        Exact Diagonalization start",16,Bold,Orange]];
		Print["e = ",e];
		Print["V = ",V];
		If[EdMode=="Superc",Print["\[CapitalDelta] = ",\[CapitalDelta]]];

		(* Build and diagonalize the AIM Hamiltonian + print timing *)
		Print["E.D. time: ",First@AbsoluteTiming[
		
			Hsectors=
				ParallelTable[
					Sum[
						impHblocks[[sectorindex,j]]*Parameters[[j]]
					,{j,1,Nparams}]
					+U*impHlocal[[sectorindex]]
				,{sectorindex,Dim@QnsSectorList}];(*list of all the hamiltonians for every sector*)

			{EgsSectorList,GsSectorList}=ParallelMap[Eigs[32],Hsectors]\[Transpose];(*find the ground state for each sector*)
			EgsSectorList=Flatten@EgsSectorList;(*correctly reshape the list *)
			GsSectorList=Replace[GsSectorList,{x_List}:>x,{0,-2}](*correctly reshape the list*)

		]," sec.\n"];

		Print["Computing ground state and Green functions..."];

		(* Compute the ground state *)
		Egs=Min@EgsSectorList;(*ground state energy (lowest of all the sectors)*)
		GsSectorIndex=
			Flatten@Position[EgsSectorList,
				_?((#>Egs-DegeneracyThreshold&&#<Egs+DegeneracyThreshold)&)
			];(*sector index where the lowest energy is obtained: if this list contains more than 1 element, there is a degeneracy*)
		DegeneracyWarning=If[Dim@GsSectorIndex>1,True,(*else*)False];(*is True if the ground state is degenerate, False otherwise*)


		GFTime=First@AbsoluteTiming[(*time for computing Green Functions and Self energies*)
		
			qnsstring = Which[
						EdMode == "Normal", ";  {n,nup} = ",
						EdMode == "Superc","; sz = "
						]; (* just a stupid output string *)
			If[!DegeneracyWarning,(*if there is NO degeneracy*)
				GsSectorIndex = First@GsSectorIndex;(*extract the number from the list*)
				GsQns = QnsSectorList[[GsSectorIndex]];(*quantum numbers {n,nup} or sz of the sector with minimal energy*)
				Gs = GsSectorList[[GsSectorIndex]];(*compute the ground state and flatten properly*)
				Print["        Ground state info:"];(*print relevant information on the ground state*)
				Print["Egs = ", Egs, qnsstring, GsQns];
				ImpurityGF=ImpurityGreenFunction[L,Nf,Egs,Gs,GsSectorIndex,QnsSectorList,Hsectors,EdMode,i\[Omega]](*impurity Green function*)
			];

			If[DegeneracyWarning,(*if there is degeneracy*)
				SetSharedVariable[ImpurityGF];(*make the impurity green function a shared variable between running kernels*)
				ImpurityGF=ConstantArray[0,{NMatsubara,2,2}];(*initialize impurity green function*)
				ParallelDo[(*parallel loop over the elements of GsSectorList*)
					GsQns=QnsSectorList[[gssectorindex]];(*quantum numbers {n,nup} of the sector with minimal energy*)
					Gs=GsSectorList[[gssectorindex]];(*compute the ground state and flatten properly*)
					Print["        Ground state info:\n",(*print relevant information on the ground state*)
					"Egs = ", Egs, qnsstring, GsQns];
					ImpurityGF+=ImpurityGreenFunction[L,Nf,Egs,Gs,gssectorindex,QnsSectorList,Hsectors,EdMode,i\[Omega]]
				,{gssectorindex,GsSectorIndex},DistributedContexts->Automatic];
			ImpurityGF=ImpurityGF/(Dim@GsSectorIndex)(*divide the green function by the number of degenerate states.*)
			];


			(*Self energy of the previous iteration*)
			\[CapitalSigma]old=If[DMFTiterator==1,ConstantArray[0,{NMatsubara,2,2}],\[CapitalSigma]];
			(* Getting the inverse non interacting impurity Green function *)
			InverseGF0=WeissField[Nbath,Parameters,EdMode,i\[Omega]];
			(* Getting the inverse interacting impurity Green function *)
			InverseGF=Inverse/@ImpurityGF;
			(* Getting the Self Energy *)
			\[CapitalSigma]=InverseGF0-InverseGF;
			(* Getting the lattice local Green function *)
			LocalGF=LocalGreenFunction["Bethe",\[CapitalSigma],EdMode,i\[Omega]];
		
		];

		Print["Green functions computed. Total time: ",GFTime," sec."];(*show evaluation time*)

		(* Print *)
		Print[Style["        Exact Diagonalization completed",16,Bold,Orange]];
		Print["----------------------------------------------------------------------------------------"];
		Print[Style["           Self Consistency start ",16,Bold,Magenta]];

		(* Self Consistency Equation *)
		Print["S.C. time: ",First@AbsoluteTiming[
		
			Which[
				EdMode == "Normal",
				(*update bath parameters minimizing the distance between hybridizations*)
				{e,V} = SelfConsistencyBethe[Nbath, LocalGF, LFit, Mixing, Parameters, EdMode, i\[Omega]];
				Parameters = Flatten@{e,V},
			(* ---------------------------------------------------------------------- *)
				EdMode == "Superc",
				{e,V,\[CapitalDelta]} = SelfConsistencyBethe[Nbath, LocalGF, LFit, Mixing, Parameters, EdMode, i\[Omega]];	
				Parameters = Flatten@{e,V,\[CapitalDelta]}
			];
		
		]," sec."];

		(* Compute error and check convergence *)
		(* function used to evaluate the error: in this case diagonal Self energy. *)
		X = \[CapitalSigma][[All,1,1]];
		Xold = \[CapitalSigma]old[[All,1,1]];
		(* compute the error *)
		error = DMFTError[X,Xold];
		(* store error *)
		AppendTo[ErrorList,{DMFTiterator,error}];

		(*Print*)
		Print[Style["           Self Consistency completed",16,Bold,Magenta]];
		Print["----------------------------------------------------------------------------------------"];
		Print["----------------------------------------------------------------------------------------"];
		Print["----------------------------------------------------------------------------------------"];

		(*Exit DMFT Loop if convengerce is reached*)
		If[DMFTiterator>DMFTMinIterations&&error<DMFTerror&&LastIteration,Break[];];
		If[error<DMFTerror,LastIteration=True,(*else*)LastIteration=False];


	,{DMFTiterator,1,DMFTMaxIterations}];


	(*           Compute observables and manage output            *)
	(* Bath parameters *)
	BathParametersFile = FolderPath<>"hamiltonian_restart_U="<>ToString[Abs@U]<>".txt";
	WriteOutput[True, BathParametersFile, "BathParameters", U,
		Which[
			EdMode == "Normal", {e,V},
			EdMode == "Superc", {e,V,\[CapitalDelta]}
		]];
	
	(* Self energy *)
	SelfEnergyFile = FolderPath<>"impSigma_U="<>ToString[Abs@U]<>".m";
	WriteOutput[StoreObservablesQ, SelfEnergyFile, "SelfEnergy",U,\[CapitalSigma]];

	(* Quasiparticle Weight *)
	ZFile = FolderPath<>"z.txt";
	z = Z[\[CapitalSigma],50,i\[Omega]];(*compute the quasiparticle weight*)
	WriteOutput[StoreObservablesQ, ZFile, "z", U, z];

	(* Order Parameter and superfluid stiffness *)
	If[EdMode == "Superc",
		\[Phi]File = FolderPath<>"phi.txt";
		\[Phi] = OrderParameter[L,Nf,GsSectorIndex,QnsSectorList,GsSectorList];(*compute order parameter from the impurity model*)
		WriteOutput[StoreObservablesQ,\[Phi]File,"\[Phi]",U,\[Phi]];
		(**)
		StiffnessFile = FolderPath<>"Ds.txt";
		Ds = SuperfluidStiffnessBethe[\[CapitalSigma], LE, i\[Omega]];
		WriteOutput[StoreObservablesQ, StiffnessFile, "Ds", U, Ds];
	];

	(* Kinetic Eenergy *) 
	KineticEnergyFile = FolderPath<>"Kinetic_Energy.txt";
	Ekin = KineticEnergyBethe[\[CapitalSigma], LE, i\[Omega]];
	WriteOutput[StoreObservablesQ, KineticEnergyFile, "Ekin", U, Ekin];

	(* Occupancies *)
	OccupancyFile = FolderPath<>"Occupancies.txt";
	occupancies = 
	Table[
		Occupancies[L, Nf, GsSectorIndex, QnsSectorList, GsSectorList, type, EdMode],
		{type,{"Double","Single","Empty"}}
	];(*compute the double, single and empty occupancies*)
	WriteOutput[StoreObservablesQ, OccupancyFile, "Occupancy", U, occupancies];

	(* Density and fluctuation *)
	DensityFile = FolderPath<>"Density.txt";
	density = ImpurityDensity[L,Nf,GsSectorIndex,QnsSectorList,GsSectorList,EdMode];
	squaredensity = SquareDensity[L,Nf,GsSectorIndex,QnsSectorList,GsSectorList,EdMode];
	fluctuationdensity = Sqrt[squaredensity-density^2];
	WriteOutput[StoreObservablesQ,DensityFile,"Density",U,{density,squaredensity,fluctuationdensity}];

	(* Spin and fluctuation *)
	SpinFile = FolderPath<>"Spin.txt";
	sz = Magnetization[L, Nf, GsSectorIndex, QnsSectorList, GsSectorList, EdMode];
	squaresz = SquareSpinZ[L, Nf, GsSectorIndex, QnsSectorList, GsSectorList, EdMode];
	fluctuationsz = Sqrt[squaresz-sz^2];
	WriteOutput[StoreObservablesQ,SpinFile,"Spin",U,{sz,squaresz,fluctuationsz}];

	(* Spectral function *)
	SpectralFunctionFile = FolderPath<>"Spectral_Function_U="<>ToString[Abs@U]<>".txt";
	spectralfunction = SpectralFunction[L,Nf,Egs,GsSectorIndex,GsSectorList,QnsSectorList,Hsectors,EdMode,\[Omega],\[Eta]];
	WriteOutput[StoreObservablesQ, SpectralFunctionFile, "SpectralFunction", U, spectralfunction];

	(* Error list *)
	ErrorFile = FolderPath<>"Error_U="<>ToString[Abs@U]<>".txt";
	WriteOutput[StoreObservablesQ, ErrorFile, "Error", U, ErrorList];

,{u,NU}]
