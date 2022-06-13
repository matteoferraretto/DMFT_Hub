(* ::Package:: *)

(* ::Subtitle:: *)
(*LOAD THE DMFT PACKAGE*)


(* ::Input::Initialization:: *)
ClearAll["Global`*"];
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
U=-0.3;(* interaction energy in units of DBethe = 1.0 *)
InitializeBathMode="Default";(*path to input file of bath parameters or "Default"*)
\[CapitalDelta]0=0.2;(* starting value of symmetry breaking field if EdMode="Superc" *)
\[Mu]=0;(*chemical potential*)

(*         INFO ON DMFT LOOPS         *)
DMFTerror=10^(-5);(*convergence threshold*)
DMFTMinIterations=1;(*minimum number of iterations*)
DMFTMaxIterations=30;(*maximum number of iterations*)
DegeneracyThreshold=10^(-6);(*two states are degenerate when |Subscript[E, 1]-Subscript[E, 2]|<DegeneracyThreshold *)
Mixing=0.3;(*mixing parameter: at iteration n set {Subscript[e, n],Subscript[V, n]}=Mixing*{Subscript[e, n-1],Subscript[V, n-1]}+(1-Mixing)*{Subscript[e, n],Subscript[V, n]} *)
LFit=1000;(*number of Matsubara frequencies used for the \[Chi]^2 fit*)
RoughIntegrationQ=False;(*set True to use a rough calculation of \[Integral]d\[Epsilon]D[\[Epsilon]]..., set False to use a more refined calculation using NIntegrate*)
LE=500;(* number of channels used to divide the interval [-DBethe,+DBethe] for the rough calculation of \[Integral]d\[Epsilon]D[\[Epsilon]]... *)
FitMode = "LocalMinimum";
FitOptions = {LFit, FitMode, Mixing};

(*      INFO ON MATSUBARA AND REAL FREQUENCIES         *)
T=0.001;(*fictitious temperature to define Matsubara frequencies*)
NMatsubara=5000;(*Total number of Matsubara frequencies*)
i\[Omega]=Table[(2*n+1)*Pi*I*T,{n,0,NMatsubara-1}];(*list of Matsubara frequencies*)
\[Omega]min=-5.;\[Omega]max=5.;(*set min and max value for the set of real frequencies*)
NReal=10000;(*number of real frequencies*)
d\[Omega]=(\[Omega]max-\[Omega]min)/NReal;(*real frequency step*)
\[Omega]=Table[\[Omega]min+n*d\[Omega],{n,0,NReal}];(*list of real frequencies*)
\[Eta]=0.025;(*small shift of the pole in the imaginary axis: this avoids singularities, but introduces an artificial broadening of the spectrum*)

(*          INFO ON OUTPUT STORAGE FILE PATH         *)
LoadHamiltonianQ = False; (*set True if you want to load the Hamiltonian from a file, False if you want to generate it on the fly and then store it on a file*)
ImpHBlocksFile = FolderPath<>"ImpHBlocks_L="<>ToString[L]<>".mx";
ImpHLocalFile = FolderPath<>"ImpHLocal_L="<>ToString[L]<>".mx";
StoreBathParametersQ = False;(*want to store bath parameters? On which file? (Very recommended for serious calculations!)*)
BathParametersFile = FolderPath<>"hamiltonian_restart_U="<>ToString[Abs@U]<>".txt";

StoreSelfEnergyQ = False;(*want to store self energy? On which file?*)
SelfEnergyFile = FolderPath<>"impSigma_U="<>ToString[Abs@U]<>".txt";
StoreSpectralFunctionQ=False;
SpectralFunctionFile = FolderPath<>"Spectral_Function_U="<>ToString[Abs@U]<>".txt";
StoreZQ = False;(*want to store z? On which file?*)
ZFile=FolderPath<>"z.txt";
Store\[Phi]Q = False;(*want to store z? On which file?*)
\[Phi]File = FolderPath<>"phi.txt";
StoreOccupancyQ = False;(*want to store occupancies? On which file?*)
OccupancyFile = FolderPath<>"Occupancies.txt";
StoreDensityQ = False;(*want to store the density and the corresponding fluctuatuion?*)
DensityFile=FolderPath<>"Density.txt";
StoreSpinQ = False;(*want to store the magnetization and spin z fluctuation?*)
SpinFile = FolderPath<>"Spin.txt";
StoreKineticEnergyQ = False;(*want to store the kinetic energy? On which file?*)
KineticEnergyFile = FolderPath<>"Kinetic_Energy.txt";
StoreStiffnessQ = False;(*want to store the superfluid stiffness?*)
StiffnessFile = FolderPath<>"Ds.txt";


(* ::Subtitle:: *)
(*DMFT LOOP*)


LastIteration=False;(*allows to do one more iteration after convergence threshold is reached*)
Converged=False;(*True if DMFT has converged, false otherwise*)
ErrorList={};(*list of DMFT errors*) 

(* GET SECTORS *)
QnsSectorList = SectorList[L, EdMode](* list of quantum numbers of all the sectors {n,nup} or sz *)
DimSectorList = DimSector[L, #, EdMode]&/@QnsSectorList(* list of dimensions of all the sectors *)
Sectors = BuildSector[L, Nf, #, EdMode]&/@QnsSectorList;(* list of all the sectors *)
SectorDispatch = Dispatch@Flatten[
		MapIndexed[{#1->#2[[1]]}&,
		QnsSectorList],1];

(* GET BATH PARAMETERS *)
Which[
	EdMode=="Normal",
	{e,V} = StartingBath[InitializeBathMode, Nbath, EdMode];
	Parameters = Flatten@{e,V},
(* ----------------------------------------------------------------- *)
	EdMode=="Superc", 
	{e,V,\[CapitalDelta]} = StartingBath[InitializeBathMode, Nbath, EdMode]; 
	\[CapitalDelta]=\[CapitalDelta]0*\[CapitalDelta]; 
	(* REMOVE!!! *)
	(*e = {-0.35, -0.0, 0.0, 0.35};
	V = {0.31, 0.13, 0.12, 0.31};
	\[CapitalDelta] = {0.02, 0.02, -0.04, 0.01};*)
	(* END REMOVE *)
	Parameters=Flatten@{e,V,\[CapitalDelta]}
];
Nparams = Length[Parameters];(* total number of parameters *)

(* GET SECTORS *)
QnsSectorList = SectorList[L,EdMode];(* list of quantum numbers of all the sectors {n,nup} or sz *)
DimSectorList = DimSector[L,#,EdMode]&/@QnsSectorList;(*list of dimensions of all the sectors*)

(* INPUT RECAP *)
Print[Style["Recap of input:",16,Bold]];
Print["Nbath: ",Nbath,". Nsectors: ", Length[QnsSectorList], ". Dim. of the largest sector: ",Max@DimSectorList];

(* GET HAMILTONIANS *)
{impHblocks,impHlocal}=GetHamiltonian[L,Nf,QnsSectorList,LoadHamiltonianQ,ImpHBlocksFile,ImpHLocalFile,EdMode];

(*                   DMFT LOOP                     *)
Do[

	(* First print *)
	Print[Style["DMFT Loop n. ",20,Bold,Red],Style[DMFTiterator,20,Bold,Red]];
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
			,{sectorindex, Length[QnsSectorList]}];(*list of all the hamiltonians for every sector*)
		
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
	DegeneracyWarning=If[Length[GsSectorIndex]>1,True,(*else*)False];(*is True if the ground state is degenerate, False otherwise*)


	GFTime=First@AbsoluteTiming[(*time for computing Green Functions and Self energies*)
		
		qnsstring = Which[
					EdMode=="Normal", ";  {n,nup} = ",
					EdMode=="Superc","; sz = "
					]; (* just a stupid output string *)
		If[!DegeneracyWarning,(*if there is NO degeneracy*)
			GsSectorIndex = First@GsSectorIndex;(*extract the number from the list*)
			GsQns = QnsSectorList[[GsSectorIndex]];(*quantum numbers {n,nup} or sz of the sector with minimal energy*)
			Gs = GsSectorList[[GsSectorIndex]];(*compute the ground state and flatten properly*)
			Print["        Ground state info:"];(*print relevant information on the ground state*)
			Print["Egs = ", Egs, qnsstring, GsQns];
			ImpurityGF=ImpurityGreenFunction[L,Nf,Egs,Gs,GsQns,Hsectors,Sectors,SectorDispatch,EdMode,i\[Omega]](*impurity Green function*)
		];

		If[DegeneracyWarning,(*if there is degeneracy*)
			SetSharedVariable[ImpurityGF];(*make the impurity green function a shared variable between running kernels*)
			ImpurityGF=ConstantArray[0,{NMatsubara,2,2}];(*initialize impurity green function*)
			ParallelDo[(*parallel loop over the elements of GsSectorList*)
				GsQns=QnsSectorList[[gssectorindex]];(*quantum numbers {n,nup} of the sector with minimal energy*)
				Gs=GsSectorList[[gssectorindex]];(*compute the ground state and flatten properly*)
				Print["        Ground state info:\n",(*print relevant information on the ground state*)
				"Egs = ", Egs, qnsstring, GsQns];
				ImpurityGF+=ImpurityGreenFunction[L,Nf,Egs,Gs,GsQns,Hsectors,Sectors,SectorDispatch,EdMode,i\[Omega]]
			,{gssectorindex,GsSectorIndex},DistributedContexts->Automatic];
		ImpurityGF=ImpurityGF/(Length[GsSectorIndex])(*divide the green function by the number of degenerate states.*)
		];


		(*Self energy of the previous iteration*)
		WeissOld = If[DMFTiterator==1,ConstantArray[0,{NMatsubara,2,2}],InverseGF0];
		(* Getting the inverse non interacting impurity Green function *)
		InverseGF0 = WeissField[Nbath,Parameters,EdMode,i\[Omega]];
		WeissNew = InverseGF0; (* <----- to be substituted by the self consistency equation *)
		(* Getting the inverse interacting impurity Green function *)
		InverseGF = Inverse/@ImpurityGF;
		(* Getting the Self Energy *)
		\[CapitalSigma] = InverseGF0 - InverseGF;
		(* Getting the lattice local Green function *)
		LocalGF = LocalGreenFunction["Bethe",\[CapitalSigma],EdMode,i\[Omega]];
		\[CapitalGamma]new = Partition[#,2]&/@({# + \[Mu],0,0,	  # - \[Mu] }\[Transpose]&/@i\[Omega]) - Inverse/@LocalGF - \[CapitalSigma];
		Print@Dimensions[InverseGF0];
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
			Parameters = Flatten@{e,V,\[CapitalDelta]};
		];
		
	]," sec."];

	(* Compute error and check convergence *)
	(* compute the error *)
	errornormal = DMFTError[WeissNew[[All,1,1]],WeissOld[[All,1,1]]];
	erroranomal = DMFTError[WeissNew[[All,1,2]],WeissOld[[All,1,2]]];
	Print[errornormal,"   ",erroranomal];
	error = (errornormal + erroranomal)/2.;

	(* store error *)
	Print[Style["Error = ",Bold],Style[error,Bold]];
	AppendTo[ErrorList,{DMFTiterator,error}];

	(*Print*)
	Print[Style["           Self Consistency completed",16,Bold,Magenta]];
	Print["----------------------------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------------------------"];
	Print["----------------------------------------------------------------------------------------"];

	(*Exit DMFT Loop if convengerce is reached*)
	If[DMFTiterator>DMFTMinIterations&&error<DMFTerror&&LastIteration,Break[];];
	If[error<DMFTerror,LastIteration=True,(*else*)LastIteration=False];


,{DMFTiterator,1,DMFTMaxIterations}]






(*           Compute observables and manage output            *)
(* Bath parameters *)
WriteOutput[StoreBathParametersQ,BathParametersFile,"BathParameters",U,
	Which[
		EdMode == "Normal", {e,V},
		EdMode == "Superc", {e,V,\[CapitalDelta]}
	]];
	
(* Self energy *)
WriteOutput[StoreSelfEnergyQ,SelfEnergyFile,"SelfEnergy",U,\[CapitalSigma]];

(* Quasiparticle Weight *)
z = Z[\[CapitalSigma],50,i\[Omega]];(*compute the quasiparticle weight*)
WriteOutput[StoreZQ,ZFile,"z",U,z];

(* Order Parameter and superfluid stiffness *)
If[EdMode == "Superc",
	\[Phi] = OrderParameter[L,Nf,GsSectorIndex,QnsSectorList,GsSectorList];(*compute order parameter from the impurity model*)
	WriteOutput[Store\[Phi]Q,\[Phi]File,"\[Phi]",U,\[Phi]];
	(**)
	Ds = SuperfluidStiffnessBethe[\[CapitalSigma], LE, i\[Omega]];
	WriteOutput[StoreStiffnessQ,StiffnessFile,"Ds",U,Ds];
];

(* Kinetic Eenergy *) 
Ekin = KineticEnergyBethe[\[CapitalSigma], LE, i\[Omega]];
WriteOutput[StoreKineticEnergyQ,KineticEnergyFile,"Ekin",U,Ekin];

(* Occupancies *)
occupancies = 
Table[
	Occupancies[L,Nf,GsSectorIndex,QnsSectorList,GsSectorList,type,EdMode],
	{type,{"Double","Single","Empty"}}
];(*compute the double, single and empty occupancies*)
WriteOutput[StoreOccupancyQ,OccupancyFile,"Occupancy",U,occupancies];

(* Density and fluctuation *)
density = ImpurityDensity[L,Nf,GsSectorIndex,QnsSectorList,GsSectorList,EdMode];
squaredensity = SquareDensity[L,Nf,GsSectorIndex,QnsSectorList,GsSectorList,EdMode];
fluctuationdensity = Sqrt[squaredensity-density^2];
WriteOutput[StoreDensityQ,DensityFile,"Density",U,{density,squaredensity,fluctuationdensity}];

(*Spin and fluctuation*)
sz = Magnetization[L, Nf, GsSectorIndex, QnsSectorList, GsSectorList, EdMode];
squaresz = SquareSpinZ[L, Nf, GsSectorIndex, QnsSectorList, GsSectorList, EdMode];
fluctuationsz = Sqrt[squaresz-sz^2];
WriteOutput[StoreSpinQ,SpinFile,"Spin",U,{sz,squaresz,fluctuationsz}]

(*Spectral function*)
spectralfunction = SpectralFunction[L, Nf, Egs, GsSectorIndex, GsSectorList, QnsSectorList, Hsectors, Sectors, SectorDispatch, EdMode, \[Omega], \[Eta]];
WriteOutput[StoreSpectralFunctionQ, SpectralFunctionFile, "SpectralFunction", U, spectralfunction]
ListPlot[
	spectralfunction,
		Joined->True,
		PlotRange->All,
		AxesLabel->{"\[Omega]","DoS"},
		Filling->Axis
]

(*Error evolution*)
ListLogPlot[
	ErrorList,
		Joined->True,
		AxesLabel->{"Iteration","Error"}
]




