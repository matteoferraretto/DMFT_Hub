(* ::Package:: *)

BeginPackage["DMFT`"]


(*       GENERAL SHORTCUTS        *)
Dim::usage = "Dim@x is equivalent to Dimensions[x][[1]]. "

Eigs::usage = "Eigs[x] returns a list with the lowest eigenvalues Egs and the corresponding eigenvector gs of the given matrix, organized as: {Egs,gs}. 
Applies Arnoldi method when the size of the matrix exceeds a certain threshold x."

StartingBath::usage = "StartingBathSuperc[InitializeBathMode, Nbath, EdMode] returns a list containing the bath parameters to start the DMFT loop.
If EdMode = ''Normal'' then the output has the form {e,V}, where e is a list of Nbath+1 elements, the frist being the impurity energy, the others being the bath energies; V is a list of Nbath elements containing the bath-impurity hybridizations.
If EdMode = ''Superc'' then the output has the form {e,V,\[CapitalDelta]}, where e and V are defined as above, and \[CapitalDelta] is the list of pairs creation (annihilation) amplitudes.
InitializeBathMode is a string with the path to the file containing the bath parameters; if it is set to ''Default'', default parameters are dropped. "

GetHamiltonian::usage = "GetHamiltonian[L_,f_,QnsSectorList_,LoadHamiltonianQ_,ImpHBlocksFile_,ImpHLocalFile_] ... "

WriteOutput::usage = "WriteOutput[condition, file, label, U, data] writes ''data'' both on screen and, if ''condition'' evaluates to True, even on ''file'' in a format which depends on ''label''. "

WriteClusterOutput::usage = "Same as WriteOutput[] but adapted to cluster output format."


(*    HILBERT SPACE    *)
BuildSector::usage = "BuildSector[L, f, qns, EdMode] creates a list of basis of the Fock subpace with given quantum numbers qns. 
If EdMode = ''Normal'', then qns={n,nup} are the particle number and the spin up particle number. If nup is set to Null, the function returns all the states with n particles and any value of nup.
If EdMode = ''Superc'', then qns=sz is the total spin-z of the fermions, i.e. nup - ndw.
The states are created in the integer represenation. L is the total number of sites, f is the number of flavours (only supports f=2 at the moment). "

SectorList::usage = "SectorList[L, EdMode] returns a list of all the quantum numbers associated to the symmetric sectors of the Hilbert space. "

DimSector::usage = " DimSector[L, qns, EdMode] gives the dimension of a given sector labeled by the quantum numbers qns. If EdMode=''Normal'', qns = {n,nup}; if EdMode=''Superc'', qns = sz = nup-ndw. L is the number of sites. "

ImpHBlocks::usage = "ImpHBlocks[L,f,QnsSectorList] gives a rank-2 tensor whose elements are sparse matrices. The first index runs over the sectors, the second runs over the 2(L-1) physical processes of the non interacting Anderson impurity model. 
Each sparse matrix represents the Hamiltonian of that specific process in that specific sector (not multiplied by the energy scale parameter). Input required: number of sites (L), number of flavors (f) and list of sectors (QnsSectorList). "

ImpHLocal::usage = "ImpHLocal[L,f,QnsSectorList] gives a rank-1 tensor whose elements are sparse matrices. The index runs over the sectors and each matrix represents the local Hubbard interaction Hamiltonian (not multiplied by U). 
Input required: number of sites (L), number of flavors (f) and list of sectors (QnsSectorList).  "


(*       GREEN FUNCTIONS        *)
Lanczos::usage = " Lanczos[H,\[Epsilon],miniter,maxiter,shift,startingvector] Performs a Lanczos algorithm on the matrix H, starting from the starting vector |startingvector>. \[IndentingNewLine]If startingvector is set to Null, a default starting vector is dropped. \[IndentingNewLine]\[Epsilon] is the convergence threshold, miniter and maxiter are the minimum and maximum number of iterations respectively.\[IndentingNewLine] shift is used to improve the efficience in the computation of the lowest eigenvalue of the Krylov matrix: if you expect many iteration to occur you can set it to a really high value, otherwise you can use shift=0  "

ImpurityGreenFunction::usage = "ImpurityGreenFunction[L, f, Egs, gs, GsSectorIndex, QnsSectorList, Hsectors, EdMode, z] computes the Green function for the interacting Anderson impurity model: <gs| \[CapitalPsi]_a \[CapitalPsi]dg_b |gs> on a set of real or complex frequencies z,
 where |gs> is the ground state of the system having energy Egs in the sector labeled by GsSectorIndex. If EdMode=''Normal'', then the spinor is \[CapitalPsi] = ( c_up, c_dw ); if EdMode = ''Superc'' the spinor is a Nambu spinor of the form \[CapitalPsi] = (c_up , cdg_dw ). 
 The list of all the sectors (QnsSectorList) and the list of the Hamiltonians for every sector (Hsectors) are also required as input. The output is a list of 2x2 matrices, each corresponding to a value of z. "

LocalGreenFunction::usage = "LocalGreenFunction[Lattice, \[CapitalSigma], EdMode, z] computes the local lattice (interacting) Green function starting from the self energy \[CapitalSigma], in the specified EdMode and in the specified Lattice (for now only Lattice=''Bethe'').
The result is a list of 2x2 matrices, each one for a specific value in the list of frequencies z (just like the input self energy).  "

SpectralFunction::usage = "SpectralFunction[L, f, Egs, GsSectorIndex, GsSectorList, QnsSectorList, Hsectors, EdMode, \[Omega], \[Eta]] returns a list of pairs {\[Omega], A(\[Omega])} corresponding to the real frequency spectral function. 
Input required: 
- general information about number of sites L and flavours f; 
- information about the ground state: ground state energy (Egs), list of ground state sector indexes (GsSectorIndex), list of all the degenerate ground states (GsSectorList) 
- information about the sectors: all the quantum numbers (QnsSectorList) and the Hamiltonian of each sector (Hsectors);
- type of impurity problem (EdMode): either ''Normal'' or ''Superc'';
- the list of real frequencies \[Omega] where you want to compute the function;
- information about the broadening \[Eta] (the Green function iscomputed slightly above the real axis at \[Omega]+i\[Eta]). "

WeissField::usage = " WeissField[Nbath, Parameters, EdMode, z] "


(*          OBSERVABLES        *)
Z::usage = "Z[\[CapitalSigma], FitCutoff, i\[Omega]] computes the quasiparticle weight z from the self-energy function \[CapitalSigma](i\[Omega]) by doing a linear fit of the low frequency values of such function. FitCutoff specifies how many points are taken for the fit.
Suggested value: FitCutoff=50.  "

ImpurityDensity::usage = "ImpurityDensity[L, f, GsSectorIndex, QnsSectorList, GsSectorList, EdMode] returns the average density of the impurity site. Besides general input (number of sites L and number of flavours f) you have to provide information about the ground 
states (all the degenerate ones). GsSectorIndex is a list with the sector indexes of the degenerate ground states, QnsSectorList is the list of all the quantum numbers of all the sectors, GsSectorList is a list of all the ground state vectors of all the sectors. "

SquareDensity::usage = "SquareDensity[L, f, GsSectorIndex, QnsSectorList, GsSectorList, EdMode] returns the average square density of the impurity site. Besides general input (number of sites L and number of flavours f) you have to provide information about the ground 
states (all the degenerate ones). GsSectorIndex is a list with the sector indexes of the degenerate ground states, QnsSectorList is the list of all the quantum numbers of all the sectors, GsSectorList is a list of all the ground state vectors of all the sectors."

Magnetization::usage = "Magnetization[L, f, GsSectorIndex, QnsSectorList, GsSectorList, EdMode] returns the average magnetization of the impurity site. Besides general input (number of sites L and number of flavours f) you have to provide information about the ground 
states (all the degenerate ones). GsSectorIndex is a list with the sector indexes of the degenerate ground states, QnsSectorList is the list of all the quantum numbers of all the sectors, GsSectorList is a list of all the ground state vectors of all the sectors."

SquareSpinZ::usage = "SquareSz[L, f, GsSectorIndex, QnsSectorList, GsSectorList, EdMode] returns the average sz^2 of the impurity site. Besides general input (number of sites L and number of flavours f) you have to provide information about the ground 
states (all the degenerate ones). GsSectorIndex is a list with the sector indexes of the degenerate ground states, QnsSectorList is the list of all the quantum numbers of all the sectors, GsSectorList is a list of all the ground state vectors of all the sectors."

Occupancies::usage = "Occupancies[L, f, GsSectorIndex, QnsSectorList, GsSectorList, type, EdMode] returns the fraction of double occupancies, single occupancies or empty sites. This can be controlled through the input string 
type=''Double'',''Single'',''Empty''. It also requires the EdMode=''Normal'' or ''Superc''. Besides general input (number of sites L and number of flavours f) you have to provide information about the ground states (all the degenerate ones). 
GsSectorIndex is a list with the sector indexes of the degenerate ground states, QnsSectorList is the list of all the quantum numbers of all the sectors, GsSectorList is a list of all the ground state vectors of all the sectors. "

KineticEnergyBethe::usage = "KineticEnergyBethe[\[CapitalSigma], LE, i\[Omega]] computes the average kinetic energy of the Hubbard model in the Bethe lattice (it includes the high frequency correction due to frequencies with n>NMatsubara).
It takes as input the self energy (\[CapitalSigma]) and the list of Matsubara frequencies (i\[Omega]). LE is the number of bins used to estimate the integral over d\[Epsilon].  " 

OrderParameter::usage = "OrderParameter[L, f, GsSectorIndex, QnsSectorList, GsSectorList] returns the superconducting order parameter \[Phi] = < c_dw c_up >. Besides general input as the number of sites (L), the number of spin states (f), the list of 
quantum numbers (QnsSectorList), the function needs information about the ground states (GsSectorList and GsSectorIndex). "

SuperfluidStiffnessBethe::usage = "SuperfluidStiffnessBethe[\[CapitalSigma], LE, i\[Omega]] returns the superfluid stiffness Ds given the self energy matrix (\[CapitalSigma]), the list of Matsubara frequencies (i\[Omega]) and the number of intervals in which the \[Epsilon] integration is divided (LE)."


(*        DMFT LOOP FUNCTIONS        *)
SelfConsistencyBethe::usage = "SelfConsistencyBethe[Nbath, LocalGF, LFit, Mixing, StartingParameters, EdMode, z] returns the list of bath parameters updated after the self-consistency procedure on the d->oo Bethe lattice.
If EdMode = ''Normal'', the output list has the form {e,V}; if EdMode = ''Superc'', the output list has the form {e,V,\[CapitalDelta]}.
 The required input is:
- general information: number of bath sites (Nbath), exact diagonalization mode (EdMode);
- the local Green function (LocalGF) and the set of Matsubara frequencies (z);
- technical details like the number of Matsubara frequencies used to perform the fit (LFit), the mixing parameter (Mixing) and the starting parameters for the minimization process (StartingParameters) in a flat list. " 

DMFTError::usage = "DMFTError[Xnew,Xold] evaluates the relative distance between the lists Xnew and Xold: Total[Abs[Xnew-Xold]]/Total[Xnew]. If Xold is set to Null it returns 1, as if Xold=0. The function prints the result on screen."


Begin["`Private`"];

(* notation facilitator: shows a number corresponding to the dimension of the first level of a list*)
Dim:=(Dimensions[#][[1]]&);

(* notation facilitator: remove empty sublists from one list*)
Clean:=Select[#,UnsameQ[#,{}]&]&;

(* compute eigenvalues *)
Eigs[MinLanczosDim_Integer]:=
	If[Length[#] >= MinLanczosDim,
		-Eigensystem[-#,1,Method->{"Arnoldi","Criteria"->"RealPart"}],
	(*else*)
		Sort[Transpose[Eigensystem[#]]][[1]]
	]&;


(* get the bath parameters to start the DMFT loop *)
StartingBath[InitializeBathMode_String, Nbath_Integer, EdMode_String]:=Module[
	{e,V,\[CapitalDelta]},
	Which[
		EdMode=="Normal",
		If[
			InitializeBathMode=="Default",
			e=Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}];
			V=Table[1.,{k,1,Nbath}],	
		(*else*)
			{e,V}=Import[InitializeBathMode,"Table"];
		];
		Return[{e,V}],
	(* ---------------------------------------------- *)
		EdMode=="Superc",
		If[
			InitializeBathMode=="Default",
			e=Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}];
			V=Table[1.,{k,1,Nbath}];
			\[CapitalDelta]=Table[1.,{k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalDelta]}=Import[InitializeBathMode,"Table"];
		];
		Return[{e,V,\[CapitalDelta]}];
	]
];


(*              HILBERT SPACE SECTORS              *)
(* integer version of basis for a single flavour *)
basis[L_Integer, m_Integer]:=FromDigits[#,2]&/@Permutations[IntegerDigits[2^m-1,2,L]];

(* integer version of the full basis for all flavours *)
BASIS[L_Integer, f_Integer, qns_]:=Flatten[Outer[{##}&,##]&@@Table[basis[L,qns[[\[Sigma]]]],{\[Sigma],1,f}],f-1];

(* builds all the states in a given sector *)
BuildSector[L_Integer, f_Integer, qns_, EdMode_String]:=Module[
	{n, Nup, nupmin, nupmax, newstates, states={}, sz, QnsList},
	Which[
		EdMode == "Normal",
		{n,Nup}=qns;
		nupmin=Max@{0,n-L};(*minimum number of spin up particles that you can put in the system*)
		nupmax=Min@{n,L};(*maximum number of up particles that you can put in the system*)
		If[Nup===Null,
			Do[
				newstates=BASIS[L,f,{nup,n-nup}];
				states=AppendTo[states,newstates];
			,{nup,nupmin,nupmax}],
		(*else*)
		states=BASIS[L,f,{Nup,n-Nup}]
		],	
	(* ----------------------------------------------------- *)	
		EdMode == "Superc",
		sz=qns;
		QnsList=Select[
					Flatten[#,1]&@
						Table[{nup,ndw},{nup,0,L},{ndw,0,L}],
							((#[[1]]-#[[2]]==sz)&&(* nup - ndw = sz *)
								Total[#]<=2L)&];(* 0 \[LessEqual] ndw \[LessEqual] L *)
		states=Flatten[BASIS[L,2,#]&/@QnsList,1]
		];
	states
];

(* list of the quantum numbers of all the sectors *)
SectorList[L_Integer, EdMode_String]:=Module[
	{QnsSectorList},
	Which[
		EdMode=="Normal",
		QnsSectorList={};
		Do[
			Do[
				AppendTo[QnsSectorList,{n,nup}];
			,{nup,Max[0,n-L],Min[n,L]}]
		,{n,0,2L}],
	(* ------------------------------------------ *)	
		EdMode=="Superc",
		QnsSectorList=Range[-L,L];
	];
	QnsSectorList
];

(* dimension of a sector with fixed L and quantum numbers qns *)
DimSector[L_Integer, qns_, EdMode_String]:=Module[
	{n, nup, sz, dim},
	Which[
		EdMode=="Normal",
		{n,nup}=qns;
		dim=Binomial[L,nup]*Binomial[L,n-nup],
	(* --------------------------------------- *)	
		EdMode=="Superc",
		sz=qns;
		dim=Total@Table[
			Binomial[L,ndw+sz]*Binomial[L,ndw]
			,{ndw,Max[-sz,0],Min[L-sz,L]}]
	];
	dim
]

(*apply cdg with a given spin to a basis state: return the integer form of the resulting basis state or 0*)
cdg[L_Integer, \[Sigma]_Integer, state:{__Integer}]:=Module[
	{binarystate},
	binarystate=IntegerDigits[#,2,L]&@state;
	If[binarystate[[\[Sigma],1]]==0,
		binarystate=ReplacePart[binarystate,{\[Sigma],1}->1];
		Return[FromDigits[#,2]&/@binarystate],
	(*else*)
		Return[0]
	];
];

(* apply c with a given spin to a basis state: return the integer form of the resulting basis state or 0 *)
c[L_Integer, \[Sigma]_Integer, state:{__Integer}]:=Module[
	{binarystate},
	binarystate=IntegerDigits[#,2,L]&@state;
	If[binarystate[[\[Sigma],1]]==1,
		binarystate=ReplacePart[binarystate,{\[Sigma],1}->0];
		Return[FromDigits[#,2]&/@binarystate],
	(*else*)
		Return[0]
	];
];



(*Counts how many fermions there are before state[[\[Sigma],i]] (excluded). If you set i=1,\[Sigma]=1 you get 0; if you set i=L+1,\[Sigma]=f you get the total number of fermions in the state*)
CountFermions[L_Integer, i_Integer, \[Sigma]_Integer, state:{__Integer}]:=
If[\[Sigma]==1,
	Sum[(IntegerDigits[#,2,L]&@(state[[\[Sigma]]]))[[k]],{k,1,i-1}],
(*else*)
	Sum[(IntegerDigits[#,2,L]&@(state[[s]]))[[k]],{s,1,\[Sigma]-1},{k,1,L}]+
	If[i!=1,
		Sum[(IntegerDigits[#,2,L]&@(state[[\[Sigma]]]))[[k]],{k,1,i-1}],
	(*else*)
	0]
];

(* sign accumulated by moving c_i\[Sigma] to the correct position when applying c_i\[Sigma]|state> or cdg_i\[Sigma]|state> *)
CSign[L_Integer, i_Integer, \[Sigma]_Integer, state:{__Integer}]:=(-1)^CountFermions[L,i,\[Sigma],state];
(* sign accumulated by moving c_i1\[Sigma]1 and c_i2\[Sigma]2 to the correct positions when applying c_i\[Sigma]1 c_i\[Sigma]2 |state> or similar pairs of operators *)
CCSign[L_Integer, i1_Integer, \[Sigma]1_Integer, i2_Integer, \[Sigma]2_Integer, state:{__Integer}]:=(-1)^(CountFermions[L,i2,\[Sigma]2,state]-CountFermions[L,i1,\[Sigma]1,state]);


(*apply cdg_\[Sigma] |gs>, where |gs> belongs to the sector (m,nup) and give the resulting vector resized to fit the dimension of the sector (n+1,nup+1) or (n+1,nup) (depending on \[Sigma])*)
ApplyCdg[L_Integer, f_Integer, \[Sigma]_Integer, gs_, qns_, Sectors_, SectorsDispatch_, EdMode_String]:=Module[
	{n,nup,sz,startingsector,finalsector,dim,sign,rules,dispatch,pos,\[Psi],thread,result},
	startingsector = Sectors[[qns/.SectorsDispatch]];
	(* build the final sector *)
	Which[
		EdMode=="Normal",
		{n,nup}=qns;
		If[(\[Sigma]==1&&nup==L)||(\[Sigma]==2&&(n-nup)==L),Return[0]];  (* trivial case *)
		finalsector=Which[
			\[Sigma]==1, Sectors[[{n+1,nup+1}/.SectorsDispatch]],
			\[Sigma]==2, Sectors[[{n+1,nup}/.SectorsDispatch]]
		],
	(* --------------------------- *)
		EdMode=="Superc",
		sz=qns;
		If[(sz==-L&&\[Sigma]==2)||(sz==L&&\[Sigma]==1),Return[0]];  (* trivial case *)
		finalsector = Which[
			\[Sigma]==1, Sectors[[(sz+1)/.SectorsDispatch]],
			\[Sigma]==2, Sectors[[(sz-1)/.SectorsDispatch]]
		];
	];
	(*evaluate the dimension*)
	dim = Length[finalsector];
	(*compute the list of signs obtained moving cdg to the correct position*)
	sign=CSign[L,1,\[Sigma],#]&/@startingsector;
	(*create a dispatch that labels all these states*)
	rules=Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1];
	dispatch=Dispatch[rules];
	(**)
	(*find the new positions where the entries of gs should go after applying cdg_spin*)
	pos=(cdg[L,\[Sigma],#]&/@startingsector)/.dispatch;
	(*remove zeroes from the newposition and cut the corresponding elements in the lists gs and newpos*)
	\[Psi]=gs*sign;
	\[Psi]=Delete[\[Psi],Position[pos,0]];
	pos=Delete[pos,Position[pos,0]];
	(*create a list of rules and define the resulting array*)
	thread=Thread[pos->#]&@\[Psi];
	result=SparseArray[thread,dim];
	result	
];

(*apply c_\[Sigma] |gs>, where |gs> belongs to the sector (m,nup) and give the resulting vector resized to fit the dimension of the sector (n-1,nup-1) or (n-1,nup) (depending on \[Sigma])*)
ApplyC[L_Integer, f_Integer, \[Sigma]_Integer, gs_, qns_, Sectors_, SectorsDispatch_, EdMode_String]:=Module[
	{n,nup,sz,startingsector,finalsector,dim,sign,rules,dispatch,pos,\[Psi],thread,result},
	startingsector = Sectors[[qns/.SectorsDispatch]];
	(* build the final sector *)
	Which[
		EdMode=="Normal",
		{n,nup}=qns;
		If[(\[Sigma]==1&&nup==0)||(\[Sigma]==2&&(n-nup)==0),Return[0]];  (* trivial case *)
		finalsector=Which[
			\[Sigma]==1, Sectors[[{n-1,nup-1}/.SectorsDispatch]],
			\[Sigma]==2, Sectors[[{n-1,nup}/.SectorsDispatch]]
		],
	(* --------------------------- *)
		EdMode=="Superc",
		sz=qns;
		If[(sz==-L&&\[Sigma]==1)||(sz==L&&\[Sigma]==2),Return[0]];
		finalsector=Which[
			\[Sigma]==1, Sectors[[(sz-1)/.SectorsDispatch]],
			\[Sigma]==2, Sectors[[(sz+1)/.SectorsDispatch]]
		];
	];
	(*evaluate the dimension*)
	dim = Length[finalsector];
	(*compute the list of signs obtained moving cdg to the correct position*)
	sign=CSign[L,1,\[Sigma],#]&/@startingsector;
	(*create a dispatch that labels all these states*)
	rules=Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1];
	dispatch=Dispatch[rules];
	(**)
	(*find the new positions where the entries of gs should go after applying cdg_spin*)
	pos=(c[L,\[Sigma],#]&/@startingsector)/.dispatch;
	(*remove zeroes from the newposition and cut the corresponding elements in the lists gs and newpos*)
	\[Psi]=gs*sign;
	\[Psi]=Delete[\[Psi],Position[pos,0]];
	pos=Delete[pos,Position[pos,0]];
	(*create a list of rules and define the resulting array*)
	thread=Thread[pos->#]&@\[Psi];
	result=SparseArray[thread,dim];
	result
];



(*     HOPPING FUNCTIONS     *)
(*gives True if hopping from j to i in sector \[Sigma] is possible, False otherwise*)
HopQ[L_Integer, i_Integer, j_Integer, \[Sigma]_Integer]:=If[(IntegerDigits[#,2,L]&@#[[\[Sigma]]])[[j]]==1&&(IntegerDigits[#,2,L]&@#[[\[Sigma]]])[[i]]==0,True,False]&;
(*selects states for which hopping j\[Rule]i in sector \[Sigma] is possible*)
HopSelect[L_Integer, i_Integer, j_Integer, \[Sigma]_Integer]:=Select[#,HopQ[L,i,j,\[Sigma]]]&;
(*put 1 in position i and 0 in position j*)
CdgiCj[i_Integer, j_Integer]:=ReplacePart[#,{i->1,j->0}]&;
(*hop from site j to site i in sector \[Sigma]*)
Hop[L_Integer, i_Integer, j_Integer, \[Sigma]_Integer]:=ReplacePart[#,\[Sigma]->FromDigits[#,2]&@CdgiCj[i,j]@IntegerDigits[#,2,L]&@#[[\[Sigma]]]]&;
(*counts how many fermions are in site "i" in state "state" considering all the flavours*)
Density[L_Integer, f_Integer, i_Integer]:=Sum[IntegerDigits[#,2,L][[k,i]],{k,1,f}]&;

(*    PAIR CREATION / ANNIHILATION FUNCTIONS     *)
(*gives True if it is possible to create a pair on site i, False otherwise*)
PairCreationQ[L_Integer, i_Integer]:=If[(IntegerDigits[#,2,L]&@#[[1]])[[i]]==0&&(IntegerDigits[#,2,L]&@#[[2]])[[i]]==0,True,False]&;
(*selects states for which hopping j\[Rule]i in flavor \[Sigma] is possible*)
PairCreationSelect[L_Integer, i_Integer]:=Select[#,PairCreationQ[L,i]]&;
(*put 1 in position i of a list*)
Cdg[i_Integer]:=ReplacePart[#,{i->1}]&;
(*create a pair of particles with spin up and down in site i and return the integer version of the states.*)
CreatePair[L_Integer, i_Integer]:=ReplacePart[#,{
	1->FromDigits[#,2]&@Cdg[i]@IntegerDigits[#,2,L]&@#[[1]],
	2->FromDigits[#,2]&@Cdg[i]@IntegerDigits[#,2,L]&@#[[2]]
	}]&;


(*      IMPURITY HAMILTONIAN           *)
(* Bath+Hybridization Hamiltonian in the case EdMode = "Normal" *)
ImpHBlocksNormal[L_Integer, f_Integer, QnsSectorList_, Sectors_, SectorsDispatch_]:=Module[
	{\[Psi],\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H={};
	Do[
		Hsector={};
		\[Psi] = Sectors[[qns/.SectorsDispatch]];
		dim = Length[\[Psi]];
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Hblock = SparseArray[{},{dim,dim}];
			Which[
				flag == "Bath",
				num = Density[L,f,j]/@\[Psi];(*local density*)
				Hblock += SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock],
			(* ----------------------- *)
				flag == "Hopping",
				Do[
					\[Psi]1 = HopSelect[L,1,j,\[Sigma]]@\[Psi];
					If[Length[\[Psi]1] == 0, Continue[];];
					\[Chi] = Hop[L,1,j,\[Sigma]]/@(\[Psi]1);
					rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
					\[CapitalSigma] = (CCSign[L,1,\[Sigma],j,\[Sigma],#]&/@\[Psi]1);
					Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
				,{\[Sigma],1,f}];
				Hblock = Hblock+Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock]
			];
		,{flag, {"Bath", "Hopping"}},{j,2,L}];
		AppendTo[H, Hsector];
	,{qns, QnsSectorList}];
	H
];

(* Bath+Hybridization Hamiltonian in the case EdMode = "Superc" *)
ImpHBlocksSuperc[L_Integer, f_Integer, QnsSectorList_, Sectors_, SectorsDispatch_]:=Module[
	{\[Psi],\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H={};
	Do[
		Hsector={};
		\[Psi] = Sectors[[sz/.SectorsDispatch]];
		dim = Length[\[Psi]];
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Hblock = SparseArray[{},{dim,dim}];
			Which[
				flag == "Bath",
				num = Density[L,f,j]/@\[Psi];(*local density*)
				Hblock += SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock],
			(* --------------------------------------------------------------- *)
				flag == "Hopping",
				Do[
					\[Psi]1 = HopSelect[L,1,j,\[Sigma]]@\[Psi];
					If[Length[\[Psi]1]==0,Continue[];];
					\[Chi] = Hop[L,1,j,\[Sigma]]/@(\[Psi]1);
					rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
					\[CapitalSigma] = (CCSign[L,1,\[Sigma],j,\[Sigma],#]&/@\[Psi]1);
					Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
				,{\[Sigma],1,f}];
				Hblock = Hblock+Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock],
			(* --------------------------------------------------------------- *)
				flag=="Superc",
				\[Psi]1=PairCreationSelect[L,j]@\[Psi];
				\[Chi]=CreatePair[L,j]/@\[Psi]1;
				rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
				\[CapitalSigma]=(CCSign[L,j,1,j,2,#]&/@\[Psi]1);
				Hblock+=SparseArray[pos->\[CapitalSigma],{dim,dim}];
				Hblock=Hblock+Hblock\[ConjugateTranspose];
				AppendTo[Hsector,Hblock];
			];
		,{flag,{"Bath","Hopping","Superc"}},{j,2,L}];
		AppendTo[H,Hsector];
	,{sz,QnsSectorList}];
	H
];

(* choose which case *)
ImpHBlocks[L_Integer, f_Integer, QnsSectorList_, Sectors_, SectorsDispatch_, EdMode_String]:=Which[
	EdMode == "Normal",	ImpHBlocksNormal[L, f, QnsSectorList, Sectors, SectorsDispatch],
	EdMode == "Superc",	ImpHBlocksSuperc[L, f, QnsSectorList, Sectors, SectorsDispatch]
];

(* Local Hamiltonian EdMode = "Normal" *)
ImpHLocalNormal[L_Integer, f_Integer, QnsSectorList_, Sectors_, SectorsDispatch_]:=Module[
	{\[Psi],H,Hsector,rules,dispatch,num},
	H={};
	Do[
		\[Psi] = Sectors[[qns/.SectorsDispatch]];
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		num = Density[L,f,1]/@\[Psi];(*local density*)
		(*interaction term*)
		Hsector = (1./2)*SparseArray@DiagonalMatrix[(num-1)^2-1/2];
		AppendTo[H, Hsector];
	,{qns, QnsSectorList}];
	H
];

(* Local Hamiltonian EdMode = "Superc" *)
ImpHLocalSuperc[L_Integer, f_Integer, QnsSectorList_, Sectors_, SectorsDispatch_]:=Module[
	{\[Psi],\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num,n,nup,e,V},
	H={};
	Do[
		\[Psi]=Sectors[[sz/.SectorsDispatch]];(* *)dim=Length[\[Psi]];
		rules=Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch=Dispatch[rules];
		num=Density[L,f,1]/@\[Psi];(*local density*)
		(*interaction term*)
		Hsector=(1./2)*SparseArray@DiagonalMatrix[(num-1)^2-1/2];
		AppendTo[H,Hsector];
	,{sz,QnsSectorList}];
	H
];

(* chose which case *)
ImpHLocal[L_Integer, f_Integer, QnsSectorList_, Sectors_, SectorsDispatch_, EdMode_String]:=Which[
	EdMode == "Normal",	ImpHLocalNormal[L, f, QnsSectorList, Sectors, SectorsDispatch],
	EdMode == "Superc",	ImpHLocalSuperc[L, f, QnsSectorList, Sectors, SectorsDispatch]
];

(* Get the Hamiltonian structure once for all *)
GetHamiltonian[L_Integer, f_Integer, QnsSectorList_, Sectors_, SectorsDispatch_, LoadHamiltonianQ_, ImpHBlocksFile_String, ImpHLocalFile_String, EdMode_String]:=Module[
	{impHblocks,impHlocal},
	If[LoadHamiltonianQ,
		Print["Getting Hamiltonians from file"];
		impHblocks=Import[ImpHBlocksFile];
		impHlocal=Import[ImpHLocalFile];
		Print["Done! Let's get started! \n"],
	(*else*)
		Print["Computing Hamiltonians..."];
		Print["Time: ",First@AbsoluteTiming[
			impHblocks=ImpHBlocks[L,f,QnsSectorList,Sectors,SectorsDispatch,EdMode];
			Export[ImpHBlocksFile,impHblocks];
			impHlocal=ImpHLocal[L,f,QnsSectorList,Sectors,SectorsDispatch,EdMode];
			Export[ImpHLocalFile,impHlocal];
		]," sec."];
		Print["Done! Let's get started! \n"];
	];
	{impHblocks,impHlocal}
];


Swap[x_,y_]:=Module[{},Return[{y,x}]];

Lanczos[H_, \[Epsilon]_Real, miniter_Integer, maxiter_Integer, shift_Integer, startingvector_]:=Module[
	{a,b,dim,a0,b1,v,w,HKrilov,E0old,E0new,nfinal},
	(* initialize array of a_n :  a[1]=a_0 , a[1+n]=a_n *)
	a=ConstantArray[0,maxiter+1];
	(*initialize array of b_n : b[n]=b_n *)
	b=ConstantArray[0,maxiter];
	(*dimension of H*)
	dim=Length[H];
	(* initialize the starting vector: if startingvector=0 drop a default v *)
	If[startingvector===Null,
		v=SparseArray[{1->1.0},{dim}],
	(*else*)
		v=startingvector
	];
	(* normalize it *)
	v=v/Norm[v];
	(* initialize a new vector w = Hv *)
	w=H . v;

	a0=(Conjugate@v) . w;
	w=w-a0*v;
	b1=Norm[w];
	a[[1]]=a0;b[[1]]=b1;
	E0old=a0;
	nfinal=maxiter;
	Do[
		If[b[[n]]<\[Epsilon]&&n>=miniter,
			nfinal=n;
			Break[];
		];
		w=w/b[[n]];(* w=Subscript[v, n] *)
		v=-b[[n]]*v;(* v=-Subscript[b, n]Subscript[v, n-1]*)
		{v,w}=Swap[v,w];
		w=w+H . v;(* w=Subscript[Hv, n]-Subscript[b, n]Subscript[v, n-1] *)
		a[[n+1]]=(Conjugate@v) . w; (*Subscript[a, n] = Subscript[v, n]Subscript[Hv, n]-Subscript[b, n]Subscript[v, n]Subscript[v, n-1]*)
		w=w-a[[n+1]]*v;
		If[n<maxiter,
			b[[n+1]]=Norm[w];
		];
		(*Build hamiltonian in the Krilov subspace*)
		HKrilov=DiagonalMatrix[
			Table[b[[k]],{k,1,n}],
		1];
		HKrilov=HKrilov+HKrilov\[ConjugateTranspose];
		HKrilov+=DiagonalMatrix[
			Table[a[[k]],{k,1,1+n}]
		];
		E0new=If[shift==0,
			Min@Eigenvalues@HKrilov,
		(*else*)	
			Eigenvalues[HKrilov-DiagonalMatrix[ConstantArray[shift,1+n]]][[1]]+shift
		];
		If[Abs[E0new-E0old]<\[Epsilon]&&n>miniter,
			nfinal=n;
			Break[];
		];
		E0old=E0new;
	,{n,1,maxiter}];
	a=Take[a,nfinal+1];
	b=Take[b,nfinal];
	Return[{E0new,a,b}]
];



(*           IMPURITY (INTERACTING) GREEN FUNCTION            *)
(* Define the Green function as a continued fraction *)
(*start defining a simple fraction*)
f[z_,{a_,b_}]:=b/(z-a);
(*fold the function f several times to get a continued fraction*)
G[z_,a_,b_]:=Fold[f,
Last@b/Last@(a-z*ConstantArray[1,(Dimensions@a)[[1]]]),Reverse@Most@Transpose@{a-z*ConstantArray[1,(Dimensions@a)[[1]]],b}];

(* compute the impurity Normal Green function both with normal and superconducting bath *)
ImpurityDiagonalGreenFunction[L_Integer, f_Integer, Egs_Real, gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_String, \[Sigma]_Integer, z_]:=Module[
	{norm,cdggs,cgs,E0,a,b,bprime,GF,H,rules,dispatch,n,nup,sz,sectorindex,\[Epsilon]=1.*10^(-13),MinLancIter=2,MaxLancIter=10^3,LancShift=0},
	Which[
		EdMode=="Normal",
		{n,nup} = GsQns;
		(*apply Cdg_\[Sigma]|gs> and C_\[Sigma]|gs> and  n = <gs| Cdg_\[Sigma] C_\[Sigma] |gs>*)
		cdggs = ApplyCdg[L,f,\[Sigma],gs/Norm[gs],{n,nup},Sectors,SectorsDispatch,"Normal"];(*add a particle*)
		cgs = ApplyC[L,f,\[Sigma],gs/Norm[gs],{n,nup},Sectors,SectorsDispatch,"Normal"];(*remove a particle*)
		norm = (Conjugate@cgs) . cgs;
	(*                                       *)	
	(*           add a particle              *)
		Which[
			\[Sigma]==1, sectorindex={n+1,nup+1}/.SectorsDispatch,
			\[Sigma]==2, sectorindex={n+1,nup}/.SectorsDispatch		
		];
		H = Hsectors[[sectorindex]];
		(*Lanczos iteration starting from the input state Cdg_\[Sigma]|gs>*)
		{E0,a,b} = Lanczos[H,\[Epsilon],MinLancIter,MaxLancIter,LancShift,cdggs/Norm[cdggs]];
		(*adapt the list b to construct the continued fraction*)
		bprime = (Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});
		(*G[z,{a0,a1,a2,a3},{1,-b1^2,-b2^2,b3^2}]*)
		GF=(1-norm)*Map[G[#+Egs,a,bprime]&,z];
	(*                                          *)
	(*           remove a particle              *)
		Which[
			\[Sigma]==1,sectorindex={n-1,nup-1}/.SectorsDispatch,
			\[Sigma]==2,sectorindex={n-1,nup}/.SectorsDispatch
		];
		H=Hsectors[[sectorindex]];
		(*Lanczos iteration starting from the input state \[Psi]*)
		{E0,a,b}=Lanczos[H,\[Epsilon],MinLancIter,MaxLancIter,LancShift,cgs/Norm[cgs]];
		(*adapt the list b to construct the continued fraction*)
		bprime=(Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});
		(*G[z,{a0,a1,a2,a3},{1,-b1^2,-b2^2,b3^2}]*)
		GF+=norm*Map[G[#-Egs,-a,bprime]&,z],
	(* ------------------------------------------------------ *)	
		EdMode=="Superc",
		sz = GsQns;
		(*apply Cdg_\[Sigma]|gs> and C_\[Sigma]|gs> and  n = <gs| Cdg_\[Sigma] C_\[Sigma] |gs>*)
		cdggs = ApplyCdg[L,f,\[Sigma],gs/Norm[gs],sz,Sectors,SectorsDispatch,"Superc"];(*add a particle*)
		cgs = ApplyC[L,f,\[Sigma],gs/Norm[gs],sz,Sectors,SectorsDispatch,"Superc"];(*remove a particle*)
		norm=(Conjugate@cgs) . cgs;
	(*                                       *)	
	(*           add a particle              *)
		Which[
			\[Sigma]==1, sectorindex=(sz+1)/.SectorsDispatch,
			\[Sigma]==2, sectorindex=(sz-1)/.SectorsDispatch
		];
		H = Hsectors[[sectorindex]];
		(*Lanczos iteration starting from the input state Cdg_\[Sigma]|gs>*)
		{E0,a,b}=Lanczos[H,\[Epsilon],MinLancIter,MaxLancIter,LancShift,cdggs/Norm[cdggs]];
		(*adapt the list b to construct the continued fraction*)
		bprime=(Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});
		(*G[z,{a0,a1,a2,a3},{1,-b1^2,-b2^2,b3^2}]*)
		GF=(1-norm)*Map[G[#+Egs,a,bprime]&,z];
	(*                                          *)
	(*           remove a particle              *)
		Which[
			\[Sigma]==1,sectorindex=(sz-1)/.SectorsDispatch,
			\[Sigma]==2,sectorindex=(sz+1)/.SectorsDispatch
		];
		H=Hsectors[[sectorindex]];
		(*Lanczos iteration starting from the input state \[Psi]*)
		{E0,a,b}=Lanczos[H,\[Epsilon],MinLancIter,MaxLancIter,LancShift,cgs/Norm[cgs]];
		(*adapt the list b to construct the continued fraction*)
		bprime=(Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});
		(*G[z,{a0,a1,a2,a3},{1,-b1^2,-b2^2,b3^2}]*)
		GF+=norm*Map[G[#-Egs,-a,bprime]&,z];
	];
	GF
];

(* compute the impurity Green function with superconducting bath *)
ImpurityGreenFunctionSuperc[L_Integer, f_Integer, Egs_Real, gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, z_]:=Module[
	{cdgup0,cdgdw0,cup0,cdw0,GFAparticle,GFAhole,GFBparticle,GFBhole,GFOparticle,GFOhole,GFO,GFPparticle,GFPhole,GFP,GFA,GFB,rules0,dispatch0,sz0,GF12,GF21},
	(* quantum number of the ground state sector *)
	sz0 = GsQns;
	(* communicate all the functions to all the kernels *)
	DistributeDefinitions[ApplyCdg,ApplyC,Lanczos,G];
	(* compute cdg_s|gs> and c_s|gs> for s = up, dw *)
	{cdgup0, cdgdw0, cup0, cdw0} = 
	With[{sz=sz0},
		Parallelize[{
			ApplyCdg[L,f,1,gs/Norm[gs],sz,Sectors,SectorsDispatch,"Superc"],
			ApplyCdg[L,f,2,gs/Norm[gs],sz,Sectors,SectorsDispatch,"Superc"],
			ApplyC[L,f,1,gs/Norm[gs],sz,Sectors,SectorsDispatch,"Superc"],
			ApplyC[L,f,2,gs/Norm[gs],sz,Sectors,SectorsDispatch,"Superc"]
		}]
	];
	(* compute all the main contributions to the Green function *)
	{GFOparticle, GFOhole, GFPparticle, GFPhole, GFAparticle, GFAhole, GFBparticle, GFBhole} = 
	With[
	{sz=sz0, cdgup=cdgup0, cdgdw=cdgdw0, cup=cup0, cdw=cdw0},
		Parallelize[{
			(*          O "Particle" contribution             *)
			Module[{Odggs,sectorindex,H,E0,a,b,bprime},
				Odggs = cdgup + cdw;(* apply Odg|gs> = (Cdg_up + C_dw)|gs> *)
				sectorindex=(sz+1)/.SectorsDispatch;(*if you create an up fermion or destroy a down fermion, you go from sz to sz+1*)
				H = Hsectors[[sectorindex]];(*Hamiltonian on that sector*)
				{E0,a,b}=Lanczos[H, 1.*10^(-16), 2, 1000, 0, Odggs/Norm[Odggs]];(* Apply Lanczos starting from Odg|gs> *)
				bprime = (Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});(* adapt the list b to construct the continued fraction *)
				((Conjugate@Odggs) . Odggs)*Map[G[#+Egs,a,bprime]&,z]
			],
			(*           O "Hole" contribution               *)
			Module[{Ogs,sectorindex,H,E0,a,b,bprime},	
				Ogs = cup + cdgdw;(* apply (C_up + Cdg_dw)|gs> = O|gs> *)
				sectorindex = (sz-1)/.SectorsDispatch;(*if you remove an up fermion or create a down fermion, you go from sz to sz-1*)
				H = Hsectors[[sectorindex]];
				{E0,a,b} = Lanczos[H, 1.*10^(-16),2,1000,0,Ogs/Norm[Ogs]];
				bprime = (Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});(*adapt the list b to construct the continued fraction*)
				((Conjugate@Ogs) . Ogs)*Map[G[#-Egs,-a,bprime]&,z]
			],
			(*         P "Particle" contribution           *)
			Module[{Pdggs,sectorindex,H,E0,a,b,bprime},
				Pdggs = cdgup + I*cdw;(* apply (Cdg_up + I*C_dw)|gs> = Pdg|gs> *)
				sectorindex = (sz+1)/.SectorsDispatch;(* if you create an up fermion or destroy a down fermion, you go from sz to sz+1 *)
				H = Hsectors[[sectorindex]];(* Hamiltonian on that sector *)
				{E0,a,b} = Lanczos[H, 1.*10^(-16), 2, 1000, 0,Pdggs/Norm[Pdggs]];(*Apply Lanczos starting from Adg|gs> *)
				bprime = (Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});(*adapt the list b to construct the continued fraction*)
				((Conjugate@Pdggs) . Pdggs)*Map[G[#+Egs,a,bprime]&,z]
			],
			(*          P "Hole" contribution             *)
			Module[{Pgs,sectorindex,H,E0,a,b,bprime},
				Pgs = cup - I*cdgdw;(*apply (C_up - I*Cdg_dw)|gs> = P|gs> *)
				sectorindex = (sz-1)/.SectorsDispatch;(*if you remove an up fermion or create a down fermion, you go from sz to sz-1*)
				H = Hsectors[[sectorindex]];
				{E0,a,b} = Lanczos[H,1.*10^(-16),2,1000,0,Pgs/Norm[Pgs]];
				bprime = (Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});(*adapt the list b to construct the continued fraction*)
				((Conjugate@Pgs) . Pgs)*Map[G[#-Egs,-a,bprime]&,z]
			],
			(*          C_up "Particle" contribution             *)
			Module[{sectorindex,H,E0,a,b,bprime},
				sectorindex = (sz+1)/.SectorsDispatch;(*if you remove an up fermion or create a down fermion, you go from sz to sz-1*)
				H = Hsectors[[sectorindex]];
				{E0,a,b} = Lanczos[H,1.*10^(-16),2,1000,0,cdgup/Norm[cdgup]];
				bprime = (Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});(*adapt the list b to construct the continued fraction*)
				((Conjugate@cdgup) . cdgup)*Map[G[#+Egs,a,bprime]&,z]
			],
			(*          C_up "Hole" contribution             *)
			Module[{sectorindex,H,E0,a,b,bprime},
				sectorindex = (sz-1)/.SectorsDispatch;(*if you remove an up fermion or create a down fermion, you go from sz to sz-1*)
				H = Hsectors[[sectorindex]];
				{E0,a,b} = Lanczos[H,1.*10^(-16),2,1000,0,cup/Norm[cup]];
				bprime = (Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});(*adapt the list b to construct the continued fraction*)
				((Conjugate@cup) . cup)*Map[G[#-Egs,-a,bprime]&,z]
			],
			(*          C_dw "Particle" contribution             *)
			Module[{sectorindex,H,E0,a,b,bprime},
				sectorindex = (sz-1)/.SectorsDispatch;(*if you remove an up fermion or create a down fermion, you go from sz to sz-1*)
				H = Hsectors[[sectorindex]];
				{E0,a,b} = Lanczos[H,1.*10^(-16),2,1000,0,cdgdw/Norm[cdgdw]];
				bprime = (Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});(*adapt the list b to construct the continued fraction*)
				((Conjugate@cdgdw) . cdgdw)*Map[G[#+Egs,a,bprime]&,-z]
			],
			(*          C_dw "Hole" contribution             *)
			Module[{sectorindex,H,E0,a,b,bprime},
				sectorindex = (sz+1)/.SectorsDispatch;(*if you remove an up fermion or create a down fermion, you go from sz to sz-1*)
				H = Hsectors[[sectorindex]];
				{E0,a,b} = Lanczos[H,1.*10^(-16),2,1000,0,cdw/Norm[cdw]];
				bprime = (Flatten@{1,-ReplacePart[b^2,-1->-(b[[-1]])^2]});(*adapt the list b to construct the continued fraction*)
				((Conjugate@cdw) . cdw)*Map[G[#-Egs,-a,bprime]&,-z]
			]
		}]
	];
	GFA = GFAparticle + GFAhole;(* GF_upup(iw) *)
	GFB = -(GFBparticle + GFBhole);(* -GF_dwdw(-iw) *)
	GFO = GFOparticle + GFOhole;
	GFP = GFPparticle + GFPhole;
	(*compute the off-diagonal part using the diagonal part and GFO, GFP*)
	GF12 = (1./2.)*(GFO-I*GFP-1.*(1-I)*(GFA+GFB));
	GF21 = (1./2.)*(GFO+I*GFP-1.*(1+I)*(GFA+GFB));
	(*return a list with NMatsubara 2x2 matrices, each of them being the G.F. at that specific frequency*)
	Partition[#,2]&/@({GFA,GF12,GF21,GFB}\[Transpose])
];

(* Return a list of 2x2 matrices representing the impurity GF in the normal spinor or Nambu spinor basis depending on EdMode *)
ImpurityGreenFunction[L_Integer, f_Integer, Egs_Real, gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_String, z_]:=Module[
	{GF,GFup,GFdw,zero},
	Which[
		EdMode=="Normal",
		GFup=ImpurityDiagonalGreenFunction[L,f,Egs,gs,GsQns,Hsectors,Sectors,SectorsDispatch,EdMode,1,z];
		zero=ConstantArray[0, Length[GFup]];
		GFdw=GFup;(* ok as long as there is spin symmetry *)
		GF=Partition[#,2]&/@({GFup,zero,zero,GFdw}\[Transpose]),
	(* --------------------------------------- *)
		EdMode=="Superc",
		GF=ImpurityGreenFunctionSuperc[L,f,Egs,gs,GsQns,Hsectors,Sectors,SectorsDispatch,z];
	];
	GF
];

(* Local Green function: only supports Bethe lattice at the moment *)
LocalGreenFunction[Lattice_String, \[CapitalSigma]_, EdMode_String, z_]:=Module[
	{Gloc,Floc,zero,LE=1000,DBethe=1.0,d\[Epsilon],LocalGF},
	Which[
		Lattice!="Bethe",
		Return@Print["Error, LocalGreenFunction only supports Bethe lattice"],
		Lattice=="Bethe",
		DoS[\[Epsilon]_]:=(2./(Pi*DBethe^2))*Sqrt[DBethe^2-\[Epsilon]^2];
	];
	d\[Epsilon]=2.*DBethe/LE;
	Which[
		EdMode=="Normal",
		Gloc=d\[Epsilon]*Sum[
			DoS[\[Epsilon]]/(z-\[Epsilon]-\[CapitalSigma][[All,1,1]])
		,{\[Epsilon],-DBethe,DBethe,d\[Epsilon]}];
		zero=ConstantArray[0, Length[Gloc]];
		LocalGF=Partition[#,2]&/@({Gloc,zero,zero,Gloc}\[Transpose]),
	(* ------------------------------------------------- *)	
		EdMode=="Superc",
		Gloc=d\[Epsilon]*Total@Table[
				DoS[\[Epsilon]]*(-z-Conjugate@\[CapitalSigma][[All,1,1]]-\[Epsilon])/(Abs[z-\[CapitalSigma][[All,1,1]]-\[Epsilon]]^2+Abs[\[CapitalSigma][[All,1,2]]]^2)
			,{\[Epsilon],-DBethe,DBethe,d\[Epsilon]}];
		Floc=-\[CapitalSigma][[All,1,2]]*d\[Epsilon]*Total@Table[
				DoS[\[Epsilon]]*(1./(Abs[z-\[CapitalSigma][[All,1,1]]-\[Epsilon]]^2+Abs[\[CapitalSigma][[All,1,2]]]^2))
			,{\[Epsilon],-DBethe,DBethe,d\[Epsilon]}];
		LocalGF=Partition[#,2]&/@({Gloc,Floc,Conjugate@Floc,-Conjugate@Gloc}\[Transpose])
	];
	LocalGF
];

(* compute the spectral function on real frequencies *)
SpectralFunction[L_Integer, f_Integer, Egs_Real, GsSectorIndex_, GsSectorList_, QnsSectorList_, Hsectors_, Sectors_, SectorDispatch_, EdMode_String, \[Omega]_, \[Eta]_Real]:=Module[
	{SpectralFunction,Gs,GsQns,d\[Omega],Nreal},
	d\[Omega] = \[Omega][[2]]-\[Omega][[1]];
	Nreal = Length[\[Omega]];
	(*initialize the spectral function*)
	SetSharedVariable[SpectralFunction];
	SpectralFunction=ConstantArray[0,Nreal];
	(*compute the spectral function for every ground state and sum up*)
	Do[(*loop over the elements of GsSectorList*)
		Gs = GsSectorList[[gssectorindex]];(*compute the ground state and flatten properly*)
		GsQns = QnsSectorList[[gssectorindex]];(*quantum numbers of the ground state sector*)
		SpectralFunction += -(1./Pi)*Tr/@Im@(ImpurityGreenFunction[L, f, Egs, Gs, GsQns, Hsectors, Sectors, SectorDispatch, EdMode, \[Omega]+I*\[Eta]]);
	,{gssectorindex,Flatten@{GsSectorIndex}}];
	(*normalize the spectral function*)
	SpectralFunction = SpectralFunction/(d\[Omega]*Total@SpectralFunction);
	{\[Omega],SpectralFunction}\[Transpose]
];


(*           HYBRIDIZATION FUNCTION             *)
(* numerical evaluation of the hybridization function *)
Hybridization[Nbath_Integer, Parameters_, EdMode_String]:=Module[
	{e,V,\[CapitalDelta],eimp=0,Hbath,Himp,Hhyb,hyb},
	Which[
		EdMode=="Normal",
		(* spinor convention: \[CapitalPsi] = ( c_{k,up} , c_{k,dw} , c_{0,up} , c_{0,dw} ) in column *)
		e=Take[Parameters,{1,Nbath}];
		V=Take[Parameters,{Nbath+1,2Nbath}];
		hyb={{0,0},{0,0}};
		Do[
			Hbath=DiagonalMatrix[{e[[k]],e[[k]]}];
			Himp=DiagonalMatrix[{eimp,eimp}];
			Hhyb=DiagonalMatrix[{V[[k]],V[[k]]}];
			hyb-=Conjugate[Hhyb] . Inverse[{{#,0},{0,#}}-Hbath] . Hhyb
		,{k,1,Nbath}],
	(* -------------------------------- *)
		EdMode=="Superc",
		(* spinor convention: \[CapitalPsi] = ( c_{k,up} , cdg_{k,dw} , c_{0,up} , cdg_{0,dw} ) in column *)
		e=Take[Parameters,{1,Nbath}];
		V=Take[Parameters,{Nbath+1,2Nbath}];
		\[CapitalDelta]=Take[Parameters,{2Nbath+1,3Nbath}];
		hyb={{0,0},{0,0}};
		Do[
			Hbath={{e[[k]],\[CapitalDelta][[k]]},{\[CapitalDelta][[k]],-e[[k]]}};
			Himp=DiagonalMatrix[{eimp,-eimp}];
			Hhyb=DiagonalMatrix[{V[[k]],V[[k]]}];
			hyb-=Conjugate[Hhyb] . Inverse[{{#,0},{0,#}}-Hbath] . Hhyb
		,{k,1,Nbath}]
	];
hyb
]&;

(* analytic evaluation of the hybridization function *)
\[CapitalGamma][Nbath_Integer, Parameters_, EdMode_String, z_]:=Module[
	{hyb,e,V,\[CapitalDelta]},
	Which[
		EdMode=="Normal",
		e=Take[Parameters,{1,Nbath}];
		V=Take[Parameters,{Nbath+1,2Nbath}];
		hyb=Total@Table[(V[[k]]^2)/(z-e[[k]]),{k,1,Nbath}],
	(* ----------------------------------------------- *)	
		EdMode=="Superc",
		e=Take[Parameters,{1,Nbath}];
		V=Take[Parameters,{Nbath+1,2Nbath}];
		\[CapitalDelta]=Take[Parameters,{2Nbath+1,3Nbath}];
		hyb=-Total@Table[
				(V[[k]]^2)*(z+e[[k]])/((-I*z)^2+e[[k]]^2+\[CapitalDelta][[k]]^2)
			,{k,1,Nbath}]
	];
	hyb
];

(* analytic evaluation of non interacting green function *)
NonInteractingGreenFunction[L_Integer, e:{__Reals}, V:{__Reals}, z_]:=1.0/(z+e[[1]]-\[CapitalGamma][L,Drop[e,1],V,z]);

(* useful functions to evaluate the non interacting Green function in the superconducting case *)
A[Nbath_Integer, Parameters_, z_]:=z-\[CapitalGamma][Nbath,Parameters,"Superc",z];
B[Nbath_Integer, Parameters_, z_]:=Module[
	{e,V,\[CapitalDelta]},
	e=Take[Parameters,{1,Nbath}];
	V=Take[Parameters,{Nbath+1,2Nbath}];
	\[CapitalDelta]=Take[Parameters,{2Nbath+1,3Nbath}];
	Total@Table[(\[CapitalDelta][[k]]*V[[k]]^2)/((-I*z)^2+e[[k]]^2+\[CapitalDelta][[k]]^2),{k,1,Nbath}]
];

(* Weiss field (inverse non interacting impurity Green function *)
WeissField[Nbath_Integer, Parameters_, EdMode_String, z_]:=
	Which[
		EdMode=="Normal",
		Partition[#,2]&/@({
			z-\[CapitalGamma][Nbath,Parameters,"Normal",z],	ConstantArray[0, Length[z]],
			ConstantArray[0, Length[z]],	z-\[CapitalGamma][Nbath,Parameters,"Normal",z]
		}\[Transpose]),
	(* ---------------------------------------------------- *)
		EdMode=="Superc",
		Partition[#,2]&/@({
			A[Nbath,Parameters,z],	-Conjugate@B[Nbath,Parameters,z],
			-B[Nbath,Parameters,z],	-A[Nbath,Parameters,-z]
		}\[Transpose]) 
	];
(*
NonInteractingNormalGreenFunction[L_,Parameters_,z_]:=A[L,Parameters,-z]/(A[L,Parameters,z]*A[L,Parameters,-z]+Abs[B[L,Parameters,z]]^2);
NonInteractingAnomalousGreenFunction[L_,Parameters_,z_]:=-(Conjugate@B[L,Parameters,z])/(A[L,Parameters,z]*A[L,Parameters,-z]+Abs[B[L,Parameters,z]]^2);
NonInteractingGreenFunction[L_,Parameters_,z_]:={
	{NonInteractingNormalGreenFunction[L,Parameters,z],NonInteractingAnomalousGreenFunction[L,Parameters,z]},
	{Conjugate@NonInteractingAnomalousGreenFunction[L,Parameters,z],-NonInteractingNormalGreenFunction[L,Parameters,-z]}
	};*)



(* SELF CONSISTENCY PROCEDURES *)
SelfConsistencyBethe[Nbath_Integer, LocalGF_, LFit_Integer, Mixing_, StartingParameters_, EdMode_String, z_]:=Module[
	{esymbols,Vsymbols,\[CapitalDelta]symbols,symbols,residue,newparameters,newe,newV,new\[CapitalDelta],lists,e,V,\[CapitalDelta],\[Chi],\[Chi]normal,\[Chi]anomalous,DBethe=1.},
	Which[
		EdMode=="Normal",
		(* extract e and V from the input list StartingParameters *)
		e=Take[StartingParameters,{1,Nbath}];
		V=Take[StartingParameters,{Nbath+1,2Nbath}];
		(* define lists of symbols {e1,e2,...} and {V1,V2,...} *)
		esymbols=Table[Symbol["e"<>ToString[i]],{i,1,Nbath}];
		Vsymbols=Table[Symbol["V"<>ToString[i]],{i,Nbath}];
		(* define a list {e1,e2,...,V1,V2,...} *)
		symbols=Flatten@{esymbols,Vsymbols};
		(*distance function: target function to minimize*)
		\[Chi][symbols_]:=
			Total@Take[#,LFit]&@(Abs[
				LocalGF[[All,1,1]]*(DBethe^2)/4-(\[CapitalGamma][Nbath,symbols,"Normal",#]&/@z)
			]^2);
		(*search a local minimum of such function*)
		{residue, newparameters}=
			FindMinimum[
				\[Chi][symbols],
				{symbols, StartingParameters}\[Transpose],
				Method->"ConjugateGradient",
				MaxIterations->700,
				AccuracyGoal->5
			];
		(*update the bath*)
		{newe,newV}={esymbols,Vsymbols}/.newparameters;
		lists={newe,newV};
		{newe,newV}=SortBy[lists\[Transpose],First]\[Transpose];
		(*if required, mix the old and the new parameters*)
		e=Mixing*e+(1-Mixing)*newe;
		V=Mixing*V+(1-Mixing)*newV;	
		Return[{e,V}],
	(* -------------------------------------------------- *)
		EdMode=="Superc",
		(* extract e, V and \[CapitalDelta] from the input list StartingParameters *)
		e=Take[StartingParameters,{1,Nbath}];
		V=Take[StartingParameters,{Nbath+1,2Nbath}];
		\[CapitalDelta]=Take[StartingParameters,{2Nbath+1,3Nbath}];
		(* define lists of symbols {e2,e3,...}, {V1,V2,...} and {\[CapitalDelta]1,\[CapitalDelta]2,...} *)
		esymbols=Table[Symbol["e"<>ToString[i]],{i,Nbath}];
		Vsymbols=Table[Symbol["V"<>ToString[i]],{i,Nbath}];
		\[CapitalDelta]symbols=Table[Symbol["\[CapitalDelta]"<>ToString[i]],{i,Nbath}];
		(* define a list {e1,e2,...,V1,V2,...,\[CapitalDelta]1,\[CapitalDelta]2,...} *)
		symbols=Flatten@{esymbols,Vsymbols,\[CapitalDelta]symbols};
		(* distance function between normal compoents (1,1) *)
		\[Chi]normal[symbols_]:=
			Total@Take[#,LFit]&@(Abs[
				LocalGF[[All,1,1]]*(DBethe^2)/4. - (\[CapitalGamma][Nbath,symbols,"Superc",#]&/@z)
			]^2);
		(* distance function between anomalous compoents (2,1) *)
		\[Chi]anomalous[symbols_]:=
			Total@Take[#,LFit]&@(Abs[
				LocalGF[[All,1,2]]*(DBethe^2)/4. + (B[Nbath,symbols,#]&/@z)
			]^2);
		(* target function to minimize (average of the two) *)
		\[Chi][symbols_]:=
			(1./2.)*(\[Chi]normal[symbols]+\[Chi]anomalous[symbols]);
		(*search the global minimum of such function*)
		{residue,newparameters}=
			FindMinimum[
				\[Chi][symbols],
				{symbols, StartingParameters}\[Transpose],
				Method->"ConjugateGradient",
				MaxIterations->1000,
				PrecisionGoal->6
			];
		(*update the bath*)
		{newe,newV,new\[CapitalDelta]}={esymbols,Vsymbols,\[CapitalDelta]symbols}/.newparameters;
		lists={newe,newV,new\[CapitalDelta]};
		{newe,newV,new\[CapitalDelta]}=SortBy[lists\[Transpose],First]\[Transpose];
		(*if required, mix the old and the new parameters*)
		e=Mixing*e+(1-Mixing)*newe;
		V=Mixing*V+(1-Mixing)*newV;
		\[CapitalDelta]=Mixing*\[CapitalDelta]+(1-Mixing)*new\[CapitalDelta];
		
		Print["Fit residue = ", residue];
		Return[{e,V,\[CapitalDelta]}]
	]
];



(*           OBSERVABLES             *)
(* Quasiparticle weight *)
\[NonBreakingSpace]Z[\[CapitalSigma]_, FitCutoff_Integer, i\[Omega]_]:=Module[
	{Selfenergy,data,a,z},
	Selfenergy = \[CapitalSigma][[All,1,1]];
	data=Take[#,FitCutoff]&@Transpose@{Im@i\[Omega],Im@Selfenergy};
	a=Fit[data,{x},x]/x;
	z=1/(1-a);
	Return[z]
];

(* Impurity density *)
ImpurityDensity[L_Integer, f_Integer, GsSectorIndex_, QnsSectorList_, GsSectorList_, EdMode_String]:=Module[
	{qns,\[Psi],gs,num,density=0},
	(*loop over all the degenerate ground states*)
	Do[
		gs = Flatten@GsSectorList[[gssectorindex]];(* ground state vector *)
		qns = QnsSectorList[[gssectorindex]];(* ground state quantum number *)
		\[Psi] = BuildSector[L,f,qns,EdMode];(* basis of the sector sz *)
		num=Density[L,f,1]/@\[Psi];(* evaluate number operator of the impurity for all the basis states *)
		density+=num . (Abs[gs]^2);(* multiply by the weights stored in gs^2 *)
	,{gssectorindex,Flatten@{GsSectorIndex}}];
	(* divide by the number of degenerate ground states to normalize *)
	density = density/(Length@Flatten@{GsSectorIndex})
];

(* Square impurity density *)
SquareDensity[L_Integer, f_Integer, GsSectorIndex_, QnsSectorList_, GsSectorList_, EdMode_String]:=Module[
	{qns,\[Psi],gs,num,squaredensity=0},
	(*loop over all the degenerate ground states*)
	Do[
		gs = Flatten@GsSectorList[[gssectorindex]];(*ground state vector*)
		qns = QnsSectorList[[gssectorindex]];(*ground state quantum number*)
		\[Psi] = BuildSector[L,f,qns,EdMode];(*basis of the sector sz*)
		num=Density[L,f,1]/@\[Psi];(*evaluate number operator of the impurity for all the basis states*)
		squaredensity+=(num^2) . (Abs[gs]^2);(*multiply by the weights stored in gs^2*)
	,{gssectorindex,Flatten@{GsSectorIndex}}];
	(*divide by the number of degenerate ground states to normalize*)
	squaredensity = squaredensity/(Length@Flatten@{GsSectorIndex})
];

(* Magnetization, i.e. < s_z > *)
Magnetization[L_Integer, f_Integer, GsSectorIndex_, QnsSectorList_, GsSectorList_, EdMode_String]:=Module[
	{ImpuritySz,qns,\[Psi],gs,sp,mag=0,k},
	(*evaluate (nup-ndw)/2 in the impurity site *)
	ImpuritySz[l_]:=(1/2)*(IntegerDigits[#,2,l][[1,1]]-IntegerDigits[#,2,l][[2,1]])&;
	(*loop over all the degenerate ground states*)
	Do[
		gs = Flatten@GsSectorList[[gssectorindex]];(*ground state vector*)
		qns = QnsSectorList[[gssectorindex]];(*ground state quantum number*)
		\[Psi]=BuildSector[L,f,qns,EdMode];(*basis of the sector sz*)
		sp=ImpuritySz[L]/@\[Psi];(*evaluate Sz operator of the impurity for all the basis states*)
		mag+=sp . (Abs[gs]^2);(*multiply by the weights stored in gs^2*)
	,{gssectorindex,Flatten@{GsSectorIndex}}];
	(*divide by the number of degenerate ground states to normalize*)
	mag=mag/(Length@Flatten@{GsSectorIndex})
]

SquareSpinZ[L_Integer, f_Integer, GsSectorIndex_, QnsSectorList_, GsSectorList_, EdMode_String]:=Module[
	{ImpuritySz,qns,\[Psi],gs,sp,impuritysquaresz=0,k},
	(*evaluate (nup-ndw)/2 in the impurity site *)
	ImpuritySz[l_]:=(1/2)*(IntegerDigits[#,2,l][[1,1]]-IntegerDigits[#,2,l][[2,1]])&;
	(*loop over all the degenerate ground states*)
	Do[
		gs = Flatten@GsSectorList[[gssectorindex]];(*ground state vector*)
		qns = QnsSectorList[[gssectorindex]];(*ground state quantum number*)
		\[Psi] = BuildSector[L,f,qns,EdMode];(*basis of the sector sz*)
		sp = ImpuritySz[L]/@\[Psi];(*evaluate Sz operator of the impurity for all the basis states*)
		impuritysquaresz+=(sp^2) . (Abs[gs]^2);(*multiply by the weights stored in gs^2*)
	,{gssectorindex,Flatten@{GsSectorIndex}}];
	(*divide by the number of degenerate ground states to normalize*)
	impuritysquaresz=impuritysquaresz/(Length@Flatten@{GsSectorIndex})
]

(* Empty, single and double occupancies *)
Occupancies[L_Integer, f_Integer, GsSectorIndex_, QnsSectorList_, GsSectorList_, type_String, EdMode_String]:=Module[
	{qns,\[Psi],gs,num,density=0,k},
	k=Which[
		type=="Double",2,
		type=="Single",1,
		type=="Empty",0
	];
	(*loop over all the degenerate ground states*)
	Do[
		gs = Flatten@GsSectorList[[gssectorindex]];(*ground state vector*)
		qns = QnsSectorList[[gssectorindex]];(*ground state quantum numbers*)
		\[Psi] = BuildSector[L,f,qns,EdMode];(*basis of the sector (n,nup)*)
		num = If[#==k,1,0]&/@(Density[L,f,1]/@\[Psi]);(*evaluate k-occupancy operator of the impurity for all the basis states*)
		density += num . (gs^2);(*multiply by the weights stored in gs^2*)
	,{gssectorindex,Flatten@{GsSectorIndex}}];
	(* divide by the number of degenerate ground states to normalize *)
	density = density/(Length@Flatten@{GsSectorIndex})
];

(* Kinetic energy, i.e. < Subscript[H, non interacting] > *)
KineticEnergyBethe[\[CapitalSigma]_, LE_Integer, i\[Omega]_]:=Module[
	{Ekin,T,d\[Epsilon],Dos,Glattice,\[CapitalSigma]0,DBethe=1.0},
	(* initialize variables *)
	Ekin = 0;
	T = (i\[Omega][[2]]-i\[Omega][[1]])/(2*Pi*I);
	d\[Epsilon] = 2.*DBethe/LE;
	DoS[\[Epsilon]_]:=(2/(Pi*DBethe^2))*Sqrt[DBethe^2-\[Epsilon]^2];
	(* compute lattice Green function *)
	Glattice[\[Epsilon]_]:=(-i\[Omega]-Conjugate@\[CapitalSigma][[All,1,1]]-\[Epsilon])/(Abs[i\[Omega]-\[CapitalSigma][[All,1,1]]-\[Epsilon]]^2+Abs[\[CapitalSigma][[All,1,2]]]^2);
	(* Self energy at the last Matsubara frequency *)
	\[CapitalSigma]0 = Last@\[CapitalSigma][[All,1,1]];
	(* Non vanishing terms of the first row of Eq. 4.12 - KineticEnergy.PDF *)
	Ekin+=2*T*d\[Epsilon]*Total@Table[
		DoS[\[Epsilon]]*\[Epsilon]*Total[2*Re@Glattice[\[Epsilon]]-2*(\[Epsilon]+\[CapitalSigma]0)/(i\[Omega]^2)]
	,{\[Epsilon],-DBethe,DBethe,d\[Epsilon]}];(*  Notice that \!\(
\*SubscriptBox[\(\[Sum]\), \(n\)]\(G\((k, 
\*SubscriptBox[\(i\[Omega]\), \(n\)])\)\)\)=2\!\(
\*SubscriptBox[\(\[Sum]\), \(n \[GreaterEqual] 0\)]\(Re[G\((k, 
\*SubscriptBox[\(i\[Omega]\), \(n\)])\)]\)\)  *)
	(*Non vanishing terms of the second row of Eq. 4.12 - KineticEnergy.PDF*)
	Ekin+=(1./(2.*T))*d\[Epsilon]*Total@Table[
		DoS[\[Epsilon]]*\[Epsilon]*(-\[Epsilon]-\[CapitalSigma]0)
	,{\[Epsilon],-DBethe,DBethe,d\[Epsilon]}];
	Re@Ekin
]

(* Superconducting order parameter, i.e. < c_dw c_up > *)
OrderParameter[L_Integer, f_Integer, GsSectorIndex_, QnsSectorList_, GsSectorList_]:=Module[
	{sz,\[Psi],dim,rules,dispatch,\[Psi]1,\[Chi],\[CapitalSigma],rows,cols,pos,gs,mask,\[Phi]},
	\[Phi]=0;(* initialize \[Phi] *)
	(*loop over all the degenerate ground states*)
	Do[
		gs=Flatten@GsSectorList[[gssectorindex]];(*ground state vector*)
		sz=QnsSectorList[[gssectorindex]];(*ground state quantum number sz *)
		\[Psi]=BuildSector[L,f,sz,"Superc"];(* basis of the sector sz *)
		dim=Length[\[Psi]];(* dimension of the sector *)
		rules=Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch=Dispatch[rules];
		\[Psi]1=PairCreationSelect[L,1]@\[Psi];(* select those states where pair creation in the impurity site is possible *)
		\[Chi]=CreatePair[L,1]/@\[Psi]1;(* apply the pair creation operator to that list *)
		\[CapitalSigma]=(CCSign[L,1,1,1,2,#]&/@\[Psi]1);
		rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
		mask=SparseArray[pos->\[CapitalSigma],{dim,dim}];
		\[Phi]+=(Conjugate@gs) . (mask . gs);
	,{gssectorindex,Flatten@{GsSectorIndex}}];
	\[Phi]=\[Phi]/(Length@Flatten@{GsSectorIndex});
	-\[Phi]
];

(* Superfluid Stiffness *)
SuperfluidStiffnessBethe[\[CapitalSigma]_,LE_Integer,i\[Omega]_]:=Module[
	{DoS,Flattice,Ds,V,T,d\[Epsilon],DBethe=1.},
	(* initialize parameters *)
	T=(i\[Omega][[2]]-i\[Omega][[1]])/(2*Pi*I);
	d\[Epsilon]=2.*DBethe/LE;
	DoS[\[Epsilon]_]:=(2./(Pi*DBethe^2))*Sqrt[DBethe^2-\[Epsilon]^2];(* Bethe lattice DoS *)
	V[\[Epsilon]_]:=(DBethe^2-\[Epsilon]^2)/3.;(* vertex function for the Bethe lattice *)
	Flattice[\[Epsilon]_]:=-\[CapitalSigma][[All,1,2]]/(Abs[i\[Omega]-\[CapitalSigma][[All,1,1]]-\[Epsilon]]^2+Abs[\[CapitalSigma][[All,1,2]]]^2);(* anomalous lattice Green function *)
	(*compute the stiffness*)
	Ds=4.*T*d\[Epsilon]*Total@Table[
		DoS[\[Epsilon]]*V[\[Epsilon]]*Total[Abs[Flattice[\[Epsilon]]]^2]
	,{\[Epsilon],-DBethe,DBethe,d\[Epsilon]}];
Re@Ds
]



(* DMFT error *)
DMFTError[Xnew_, Xold_, EdMode_String]:=Module[
	{error = 0.0},
	If[
	Xold === Null, error=1.0,
	(*else*)
	Which[
		EdMode == "Normal",
		error = Total[Abs[
					Xnew[[All, 1, 1]] - Xold[[All, 1, 1]]
				]]/Max[
					Total[Abs[Xnew[[All, 1, 1]]]], Total[Abs[Xold[[All, 1, 1]]]]
			],
		(* ------------------------------- *)
		EdMode == "Superc",
		error = Mean[{
					Total[Abs[
						Xnew[[All, 1, 1]] - Xold[[All, 1, 1]]
					]]/Max[
						Total[Abs[Xnew[[All, 1, 1]]]], Total[Abs[Xold[[All, 1, 1]]]]		
				],
				Total[Abs[
						Xnew[[All, 1, 2]] - Xold[[All, 1, 2]]
					]]/Max[
						Total[Abs[Xnew[[All, 1, 2]]]], Total[Abs[Xold[[All, 1, 2]]]]		
				]
			}]
	];
	];
	error
];

(* general template for output storage *)
WriteOutput[condition_, file_String, label_String, U_Real, data_]:=Module[
	{fout},
	Which[
		label=="BathParameters",
		If[condition,
			Export[file,data,"Table"];
			Print["Converged bath parameters stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label=="SelfEnergy",
		If[condition,
			Export[file,data,"Table"];
			Print["Self Energy stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label=="Error",
		If[condition,
			Export[file,data,"Table"];
			Print["Error list stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label=="SpectralFunction",
		If[condition,
			Export[file,data,"Table"];
			Print["Spectral Function stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label=="z",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data,"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["z = ",data,"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="\[Phi]",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data,"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["\[Phi] = ",data,"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Ekin",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data,"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["\!\(\*SubscriptBox[\(E\), \(kin\)]\) = ",data,"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Ds",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data,"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["\!\(\*SubscriptBox[\(D\), \(s\)]\) = ",data,"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Occupancy",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]]," ",data[[2]]+2.*data[[3]],"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["U","\t\t\t","Double Occ.","\t\t\t","Single Occ.","\t\t\t","Empty Occ.","\t\t\t","Density"];
		Print[U,"\t\t\t",data[[1]],"\t\t\t ",data[[2]],"\t\t\t\t",data[[3]],"\t\t\t\t",data[[2]]+2.*data[[3]],"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Density",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]],"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["U","\t\t\t\t","<n>","\t\t\t\t","<\!\(\*SuperscriptBox[\(n\), \(2\)]\)>","\t\t\t\t","\[CapitalDelta]n"];
		Print[U,"\t\t\t",data[[1]],"\t\t\t ",data[[2]],"\t\t\t\t",data[[3]],"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Spin",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]],"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["U","\t\t\t\t","<\!\(\*SubscriptBox[\(s\), \(z\)]\)>","\t\t\t\t","<\!\(\*SuperscriptBox[SubscriptBox[\(s\), \(z\)], \(2\)]\)>","\t\t\t","\!\(\*SubscriptBox[\(\[CapitalDelta]s\), \(z\)]\)"];
		Print[U,"\t\t\t",data[[1]],"\t\t\t ",data[[2]],"\t\t\t\t",data[[3]],"\n"](* print data on screen *)
	]
];

(* general template for output storage *)
WriteClusterOutput[condition_, file_String, label_String, U_Real, data_]:=Module[
	{fout},
	Which[
		label=="BathParameters",
		If[condition,
			Export[file,data,"Table"];
			Print["Converged bath parameters stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label=="SelfEnergy",
		If[condition,
			Export[file,data,"Table"];
			Print["Self Energy stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label=="Error",
		If[condition,
			Export[file,data,"Table"];
			Print["Error list stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label=="SpectralFunction",
		If[condition,
			Export[file,data,"Table"];
			Print["Spectral Function stored on file.\n"]
		],
	(* ---------------------------------------- *)
		label=="z",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data,"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["z = ",data,"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="\[Phi]",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data,"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["phi = ",data,"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Ekin",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data,"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["Ekin = ",data,"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Ds",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data,"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["Ds = ",data,"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Occupancy",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]]," ",data[[2]]+2.*data[[3]],"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["U","\t\t","Double Occ.","\t\t","Single Occ.","\t\t","Empty Occ.","\t\t","Density"];
		Print[U,"\t\t",data[[1]],"\t\t ",data[[2]],"\t\t\t",data[[3]],"\t\t",data[[2]]+2.*data[[3]],"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Density",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]],"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["U","\t\t","<n>","\t\t","<n2>","\t\t","Delta n"];
		Print[U,"\t\t",data[[1]],"\t\t ",data[[2]],"\t\t",data[[3]],"\n"],(* print data on screen *)
	(* ---------------------------------------- *)
		label=="Spin",
		If[condition,
			fout=OpenAppend[file];(* open output stream on file *)
			WriteString[fout,U," ",DecimalForm@data[[1]]," ",data[[2]]," ",data[[3]],"\n"]; (*write the data string*)
			Close[fout](* close output stream *)
		];
		Print["U","\t\t","<sz>","\t\t","<sz2>","\t\t","Delta sz"];
		Print[U,"\t\t",data[[1]],"\t\t ",data[[2]],"\t\t",data[[3]],"\n"](* print data on screen *)
	]
];


PlotState[L_Integer, Norb_Integer, state_]:=Module[
	{
	uparrow={Arrowheads[Large],Arrow[{{0,0},{0,.25}}]},
	downarrow={Arrowheads[Large],Arrow[{{.25,.25},{.25,0}}]},
	colors=Table[Hue[rgb],{rgb,1/Norb,1,1/Norb}],
	gup,gdw,gdocc,gempty,grid,
	binarystate
	},
	binarystate=IntegerDigits[#,2,L]&@state;
	If[Norb==1,binarystate={binarystate}];
	gup[color_]:=Graphics[{{color,uparrow},Circle[{.25,.125},.02]}];
	gdw[color_]:=Graphics[{Circle[{0,.125},.02],{color,downarrow}}];
	gdocc[color_]:=Graphics[{{color,uparrow},{color,downarrow}}];
	gempty[color_]:=Graphics[{Circle[{0,.125},.02],Circle[{.25,.125},.02]}];
	grid=Table[
		Which[
			binarystate[[iorb,1,j]]==1&&binarystate[[iorb,2,j]]==1,
			gdocc[colors[[iorb]]],
		(* ------------- *)
			binarystate[[iorb,1,j]]==1&&binarystate[[iorb,2,j]]==0,
			gup[colors[[iorb]]],
		(* ------------- *)
			binarystate[[iorb,1,j]]==0&&binarystate[[iorb,2,j]]==1,
			gdw[colors[[iorb]]],
		(* ------------ *)
			binarystate[[iorb,1,j]]==0&&binarystate[[iorb,2,j]]==0,
			gempty[colors[[iorb]]]
		]
		,{iorb,1,Norb},{j,1,L}];
	GraphicsGrid[grid,Frame->All,AspectRatio->1/3]
];

End[];

EndPackage[];
