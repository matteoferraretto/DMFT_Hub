(* ::Package:: *)

BeginPackage["DMFT`"]


(* General facilities *)
Eigs::usage = "Eigs[H] returns the lowest eigenvalue(s) and the corresponding eigenvector(s) of H in the form {values, vectors}. There are several optional arguments: '
'Temperature'', ''MinLanczosDim'', ''DegeneracyThreshold'', ''BoltzmannThreshold''. 
If ''Temperature'' -> 0 (default), then the function returns only the lowest eigenstate(s) with the respective degeneracies. Two states are considered degenerate when the difference of 
their energies is below the desired ''DegeneracyThreshold'' (\!\(\*SuperscriptBox[\(10\), \(-9\)]\) is the default).
If ''Temperature'' is > 0, then all the lowest eigenstates with a Boltzmann weight above ''BoltzmannThreshold'' are returned (\!\(\*SuperscriptBox[\(10\), \(-9\)]\) is the default). 
The option ''MinLanczosDim'' establishes a threshold on the matrix dimension, above which the Lanczos method is adopted and below which full diagonalization is performed (32 by default)."

StartingBath::usage = "StartingBath[L, f, Norb, InitializeBathMode, EdMode] returns a list containing the bath parameters to start the DMFT loop.
If EdMode = ''Normal'' then the output has the form {e,V}, where e and V are lists of Norb x (L-1) elements, representing the bath energies and the bath-impurity hybridizations.
If EdMode = ''Superc'' then the output has the form {e,V,\[CapitalDelta]}, where e and V are defined as above, and \[CapitalDelta] is the Norb x Nbath dimensional list of pairs creation (annihilation) amplitudes.
If EdMode = ''InterorbSuperc'' then the output has the form {e,V,\[CapitalDelta],\[CapitalXi]}, where e, V, \[CapitalDelta] are as above, and \[CapitalXi] is the Nbath - dimensional list of interorbital pairs creation (annihilation) amplitudes
InitializeBathMode is a string with the path to the file containing the bath parameters; if it is set to ''Default'', default parameters are dropped."

TakeIndependentParameters::usage = "TakeIndependentParameters[L, f, Norb, \[Sigma], orb, BathParameters, EdMode] returns a flat list of independent bath parameters depending on EdMode.
The word ''independent'' means that in case there is some symmetry, for example orbital and spin symmetry, we just extract a representative subset of the bath parameters with labels 
\[Sigma] and orb, namely e_\[Sigma],orb ; V_\[Sigma],orb, etc. This is useful in two cases: to compute the numerical Weiss field from the symbolic expression, and to perform the self consistency minimization
with the smallest possible number of variables. "

ReshapeBathParameters::usage = "ReshapeBathParameters[L, f, Norb, IndependentParameters, EdMode] takes a list of independent bath parameters and returns a reshaped version of it which is 
consistent with the conventionally chosen shape. "


(* sector defining functions *)
SectorList::usage = "SectorList[L, f, Norb, EdMode] returns a list of quantum numbers labeling the sectors which divide the full Hilbert space. 
If EdMode = ''Normal'' each element in the list has the form {n_orb1_spin1, n_orb1_spin2, ..., n_orb2_spin1, n_orb2_spin2, ...}
If EdMode = ''InterorbNormal'' each element in the list has the form {n_spin1, n_spin2, ...} (total number of particles with given spin)
If EdMode = ''Superc'' each element in the list has the form {sz_orb1, sz_orb2, ...} where sz_orb = n_orb_up-n_orb_dw (it only makes sense for spin 1/2 particles so far).
If EdMode = ''InterorbSuperc'' each element in the list has the form sz_tot where sz_tot = \!\(\*SubscriptBox[\(\[Sum]\), \(orb\)]\)(n_orb_up - n_orb_dw) (it only makes sense for spin 1/2 particles so far)."

DimSector::usage = "DimSector[L, f, Norb, qns, EdMode] returns the theoretical value of the dimension of a sector labeled by the quantum number(s) qns."

BuildSector::usage = "BuildSector[L, f, Norb, qns, EdMode] creates a list of basis states of the Fock subpace with given quantum numbers qns. 
If EdMode = ''Normal'', then qns = {n_orb1_spin1, n_orb1_spin2, ..., n_orb2_spin1, n_orb2_spin2, ...}.
If EdMode = ''InterorbNormal'', then qns = {n_spin1, n_spin2, ...} where n_spin is the total number of particles with a given spin.
If EdMode = ''Superc'', then qns = {sz_orb1, sz_orb2, ...} are the total spin-z in the given orbital, i.e. sz_orb = n_orb_up - n_orb_dw.
If EdMode = ''InterorbSuperc'', then qns = sz is the total spin-z.
The states are created in the integer represenation. L is the total number of sites, f is the number of flavours and Norb is the number of orbitals. 
The superconductive EdModes only support f=2 at the moment. "

DrawState::usage = "DrawState[L, f, Norb] draws a graphic representation of a Fock state that can be manipulated. Each box can be filled either with 0 (no particles in that slot) or 1 (a particle in that slot). "


(* Impurity Hamiltonian defining functions *)
HNonlocal::usage = "HNonlocal[L, f, Norb, Sectors, EdMode]"

HNonlocalInfo::usage = "HNonlocalInfo[L, f, Norb, EdMode] prints useful info about the order of non local Hamiltonian blocks."

HLocal::usage = "HLocal[L, f, Norb, Sectors, EdMode] "

GetHamiltonian::usage = "GetHamiltonian[L_, f_, Norb_, Sectors_, LoadHamiltonianQ_, HnonlocFile_, HlocFile_, EdMode_]"

HImp::usage = "HImp[Norb_, HnonlocBlocks_, HlocBlocks_, BathParameters_, InteractionParameters_, EdMode_]"

cdg::usage = "."
c::usage = "."
CSign::usage = "."
CreateParticleSelect::usage = "."
DestroyParticleSelect::usage = "."
GreenFunctionED::usage = "GreenFunctionED[L, f, Norb, {i,j}, \[Sigma], orb, Sectors, QnsSectorList, eigs, T, zlist, EdMode]"
GreenFunctionImpurity::usage = "GreenFunctionImpurity[L_, f_, Norb_, \[Sigma]_, orb_, Egs_, gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_]"
GreenFunctionImpurityNambu::usage = "GreenFunctionImpurityNambu[L_, f_, Norb_, orb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_]"
LocalGreenFunction::usage = "."
InverseGreenFunction::usage = "InverseGreenFunction[L, f, Norb, \[Sigma], orb, Egs, gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist] evaluates numerically the inverse Green function 
of the fully interacting impurity problem. This function returns a list of either numbers when EdMode = Normal, or matrices in the suitable Nambu basis depending on EdMode. The input 
spin and orbital indexes are completely ignored when the Nambu basis mixes up spin / orbital degrees of freedom. "
GreenFunction0::usage = "GreenFunction0[L_, f_,\[Mu]_, symbols_, z_, EdMode_]"
WeissField::usage = "WeissField[L, f, \[Mu], symbols, z, EdMode] returns a symbolic expression for the Weiss field of the problem. The input variable ''symbols'' should be a list of symbols
with the minimal number of bath parameters needed to build the Weiss field for the given problem. For example, if L=2, f=2, Norb=2, EdMode = ''Normal'' and the problem has orbital and spin 
symmetry, then symbols = {e1, e2, V1, V2} [do not provide the whole list of 2*L*f*Norb bath parameters, but only the 2L representatives]. "


(* Fermionic ladder Hamiltonian defining functions *)
Hnonint::usage = "Hnonint[L, f, Norb, Sectors, EdMode]. Optional arguments: RealPBC -> True (default)/ False; SyntheticPBC -> True / False (default), RealPhase -> {0}, SyntheticPhase -> 0  "


(* Self consistency *)
SelfConsistency::usage = "SelfConsistency[DBethe_, \[Mu]_, Weiss_, symbols_, StartingParameters_, LocalG_, zlist_, EdMode_]"


Begin["`Private`"];


(* compute eigenstates *)
Eigs[H_, OptionsPattern[]] := Module[
	{values, vectors, eigs, 
	dim = Length[H], 
	T = OptionValue["Temperature"], 
	\[Epsilon]deg = OptionValue["DegeneracyThreshold"], 
	\[Epsilon]temp = - OptionValue["Temperature"] * Log[OptionValue["BoltzmannThreshold"]]
	},
	If[dim <= OptionValue["MinLanczosDim"],
		(* if the matrix is small, just compute the full spectrum *)
		eigs = Eigensystem[H];
		eigs = SortBy[eigs\[Transpose], First]\[Transpose];
		values = eigs[[1]];
		(* return a suitable number of states depending on temperature *)
		If[T == 0,
			(* return just the ground state, taking into account possible degeneracies *)
			Return[
				Take[#,
					Count[
						values, x_/;(Abs[x - values[[1]]] < \[Epsilon]deg)
					]
				]&/@eigs
			],
		(* else, if T != 0 *)
			Return[
				Take[#,
					Count[
						values, x_/;((Abs[x - values[[1]]] < \[Epsilon]temp) || (Abs[x - values[[1]]] < \[Epsilon]deg))
					]
				]&/@eigs
			]
		],
	(* else *)
		(* if the matrix is large, apply Lanczos *)
		Do[
			eigs = -Eigensystem[-H, n, Method->{"Arnoldi","Criteria"->"RealPart"}];
			values = eigs[[1]];
			(* sometimes if H has degeneracies, eigenvalues are not sorted properly. We fix this by the following code *)
			If[values[[1]] > values[[2]], values = Sort[values]; eigs = SortBy[eigs\[Transpose], First]\[Transpose];];
			(* end of fixing line *)
			(* condition for stopping the loop *)
			If[T == 0,
				(* if there are no degeneracies, break the loop *)
				If[Abs[values[[-1]] - values[[-2]]] > \[Epsilon]deg, Break[];];,
			(* else if T != 0 *)
				(* if the n-th Boltzmann weight is below threshold and there are no degeneracies, break the loop *)
				If[(Abs[values[[-1]] - values[[-2]]] > \[Epsilon]deg) && (Abs[values[[-1]] - values[[1]]] > \[Epsilon]temp), Break[];];
			];
		, {n, 2, dim}];
		Return[Drop[#,-1] &/@ eigs]
	];
];
Options[Eigs] = {"Temperature" -> 0, "MinLanczosDim" -> 32, "DegeneracyThreshold" -> 10^(-9), "BoltzmannThreshold" -> 10^(-9)};

(* Initialize starting bath *)
StartingBath[L_, f_, Norb_, InitializeBathMode_, EdMode_, OptionsPattern[]] := Module[
	{e, V, \[CapitalDelta], \[CapitalXi], Nbath},
	Nbath = L - 1;
	Which[
		EdMode == "Normal" || EdMode == "InterorbNormal",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb],	
		(*else*)
			{e,V} = Import[InitializeBathMode,"Table"];
		];
		Return[{e, V}],
(* ---------------------------------------------- *)
		EdMode == "Superc",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb];
			\[CapitalDelta] = OptionValue[\[CapitalDelta]0] *ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb],
		(*else*)
			{e,V,\[CapitalDelta]} = Import[InitializeBathMode,"Table"];
		];
		Return[{e, V, \[CapitalDelta]}],
(* ---------------------------------------------- *)
		EdMode == "InterorbSuperc",
		If[
			InitializeBathMode=="Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb];
			V = ConstantArray[
				Table[1., {k,1,Nbath}],
			f*Norb];
	        \[CapitalXi] = OptionValue[\[CapitalXi]0] * Table[1., {k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalXi]} = Import[InitializeBathMode,"Table"];
		];
	   Return[{e,V,\[CapitalXi]}],
(* ---------------------------------------------- *)
		EdMode == "FullSuperc",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			f*Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			f*Norb];
			\[CapitalDelta] = OptionValue[\[CapitalDelta]0] * ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb];
	        \[CapitalXi] = OptionValue[\[CapitalXi]0] * Table[1.,{k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalDelta],\[CapitalXi]} = Import[InitializeBathMode,"Table"];
		];
	   Return[{e,V,\[CapitalDelta],\[CapitalXi]}];
	]
];
Options[StartingBath] = {\[CapitalDelta]0 -> 1., \[CapitalXi]0 -> 1.};

EdModeInfo[EdMode_] := Which[
	EdMode == "Normal",
	Print["The sectors' quantum numbers are the number of fermions for each flavor and each orbital."],
	EdMode == "Superc",
	Print["The sectors' quantum numbers are the total spin_z operators for each orbital. The bath can exchange pairs with a reservoir, but pairs have an orbital index."],
	EdMode == "InterorbNormal",
	Print["The sectors' quantum numbers are the total number of fermions for each flavor (NOT orbital-wise)."],
	EdMode == "Raman",
	Print["The sectors' quantum numbers are the orbital-wise total number of fermions (NOT flavor-wise)"],
	EdMode == "InterorbSuperc",
	Print["The sectors' quantum numbers are the total spin_z operators (NOT orbital-wise). The bath can exchange pairs with a reservoir, but pairs are inherently interorbital."],
	EdMode == "FullSuperc",
	Print["The sectors' quantum numbers are the total spin_z operators (NOT orbital-wise). The bath can exchange pairs with a reservoir, but pairs are both intraorbital and interorbital."]
];


(*              HILBERT SPACE SECTORS              *)
(* integer version of basis for a single flavour *)
basis = Compile[{
	{L, _Integer}, {m, _Integer}
	},
	FromDigits[#,2]&/@Permutations[IntegerDigits[2^m-1,2,L]],
	RuntimeAttributes->{Listable},Parallelization->True
];

(* integer version of the full basis for all flavours *)
BASIS[L_, flavors_, qns_] := Flatten[Outer[{##}&,##]&@@Table[basis[L,qns[[\[Sigma]]]],{\[Sigma],1,flavors}],flavors-1];

(* list of all the quantum numbers that label the sectors *)
SectorList[L_, f_, Norb_, EdMode_]:=Module[
	{QnsSectorList},
	Which[
		EdMode == "Normal",
		QnsSectorList=Flatten[
		Outer[{##}&,##]&@@ConstantArray[
			Range[0,L]
		,Norb*f],Norb*f-1],
(* --------------------------------------- *)
		EdMode == "InterorbNormal",
		QnsSectorList=Flatten[
		Outer[{##}&,##]&@@ConstantArray[
			Range[0,Norb*L]
		,f],f-1],
(* --------------------------------------- *)		
		EdMode == "Raman",
		QnsSectorList = Flatten[
		Outer[{##}&,##]&@@ConstantArray[
			Range[0,f*L]
		,Norb], Min[1,Norb-1]],
(* --------------------------------------- *)
		EdMode == "Superc",
		QnsSectorList=Flatten[
		Outer[{##}&,##]&@@ConstantArray[
			Range[-L,L]
		,Norb], Min[1,Norb-1]],
(* --------------------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		QnsSectorList = Range[-Norb*L,Norb*L];
	];
	QnsSectorList
];

(* return the theoretical estimate of the dimension of a sector labeled by the quantum number(s) qns *)
DimSector[L_, f_, Norb_, qns_, EdMode_] := Module[
	{n, nup, sz, dim},
	Which[
		EdMode == "Normal",
		dim = Product[
			Binomial[L,qns[[i]]]
		,{i,1,Norb*f}],
(* --------------------------------------- *)
		EdMode == "InterorbNormal",
		dim = Product[
			Binomial[Norb*L,qns[[i]]]
		,{i,1,f}],
(* --------------------------------------- *)	
		EdMode == "Superc",
		dim = Product[
			Total@Table[
				Binomial[L,ndw+qns[[i]]]*Binomial[L,ndw]
			,{ndw,Max[-qns[[i]],0],Min[L-qns[[i]],L]}]
		,{i,1,Norb}],
(* --------------------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		dim = Total@Table[
			Binomial[Norb*L,ndw+qns]*Binomial[Norb*L,ndw]
		,{ndw,Max[-qns,0], Min[Norb*L-qns,Norb*L]}];
	];
	dim
];

(* build sector, i.e. list of all the Fock states with given quantum number(s) qns *)
BuildSector[L_, f_, Norb_, qns_, EdMode_] := Module[
	{QnsList,states},
	Which[
		EdMode == "Normal",
		states = BASIS[L,f*Norb,qns],
(* ---------------------------------------------- *)
		EdMode == "InterorbNormal",
		QnsList = SectorList[L,f,Norb,"Normal"];
		QnsList = Select[
			QnsList,
			({Sum[#[[2i-1]],{i,1,Norb}],Sum[#[[2i]],{i,1,Norb}]}== qns)&
		];
		states = Flatten[BASIS[L,f*Norb,#]&/@QnsList, 1],
(* ---------------------------------------------- *)
		EdMode == "Superc",
		QnsList = SectorList[L,f,Norb,"Normal"];
		QnsList = Select[
			QnsList,
			(Delete[#,Table[{2*i},{i,1,Norb}]]-Delete[#,Table[{2*i-1},{i,1,Norb}]] == qns)&
		];
		states = Flatten[BASIS[L,Norb*f,#]&/@QnsList,1],
(* ---------------------------------------------- *)
		EdMode == "Raman",
		QnsList = SectorList[L,f,Norb,"Normal"];
		QnsList = Select[
			QnsList,
			(Total/@Partition[#,f]== qns)&
		];
		states = Flatten[BASIS[L,Norb*f,#]&/@QnsList,1],
(* ---------------------------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		QnsList = SectorList[L,f,Norb,"Normal"];
		QnsList = Select[
			QnsList, 
	    	(Table[(-1)^(1+j),{j,1,Norb*f}] . #==qns)&
	    ];
		states = Flatten[BASIS[L,Norb*f,#]&/@QnsList,1]
	];
	states
];

(* Draw a picture of a state to help the user *)
DrawState[L_, f_, Norb_, j_, \[Sigma]_, orb_] := Module[
	{impuritycolor, color},
	impuritycolor[i_]:=If[i==1,Blue,Black];
	color[i_,index_]:=If[i==j && index==f*(orb-1)+\[Sigma], Red, Black];
	Table[
		Graphics[
			Table[
				{EdgeForm[{Thick,impuritycolor[i]}],color[i,index],Opacity[.3],Rectangle[{1.1*i,0}]},{i,L}
			]
		],
	{index, f*Norb}]
];
DrawState[L_, f_, Norb_] := Module[{},
	Print[Style["Architecture of a state",16]];
	Print[Style["Red:",Red]," (site j, spin \[Sigma], orbital orb)"];
	Print[Style["Blue edge:",Blue]," impurity"];
	Manipulate[
		DrawState[L,f,Norb,j,\[Sigma],orb],
	{j,1,L,1}, {\[Sigma],1,f,1}, {orb,1,Norb,1}]
];


(*               MANIPULATION OF STATES                *)
(* gives True if it is possible to create a particle on site (i,\[Sigma],orb), False otherwise *)
CreateParticleQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}
	},
	(IntegerDigits[#,2,L]&@state[[f*(orb-1)+\[Sigma]]])[[i]] == 0,
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True
];

(* gives True if it is possible to destroy a particle on site (i,\[Sigma],orb), False otherwise *)
DestroyParticleQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}
	},
	(IntegerDigits[#,2,L]&@state[[f*(orb-1)+\[Sigma]]])[[i]] == 1,
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True
];

(* select states for which it is possible to create a particle on site (i,\[Sigma],orb) *)
CreateParticleSelect[L_, f_, i_, \[Sigma]_, orb_, stateList_] := Module[
	{criteria = CreateParticleQ[L, f, i, \[Sigma], orb, stateList]},
	Pick[stateList, criteria]
];

(* select states for which it is possible to destroy a particle on site (i,\[Sigma],orb) *)
DestroyParticleSelect[L_, f_, i_, \[Sigma]_, orb_, stateList_] := Module[
	{criteria = DestroyParticleQ[L, f, i, \[Sigma], orb, stateList]},
	Pick[stateList, criteria]
];

(* apply cdg_orb_spin in the impurity site to a basis state: return the integer form of the resulting basis state or 0 *)
cdg = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}
	},
	MapAt[
		BitSet[#, L-i]&,
		state,
		f*(orb-1)+\[Sigma]
	], CompilationTarget->"C", RuntimeAttributes->{Listable}];

(* apply c_orb_spin in the impurity site to a basis state: return the integer form of the resulting basis state or 0 *)
c = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}
	},
	MapAt[
		BitClear[#, L-i]&,
		state,
		f*(orb-1)+\[Sigma]
	], CompilationTarget->"C", RuntimeAttributes->{Listable}];

(* Counts how many fermions there are before state[[\[Sigma],orb,i]] (excluded). If you set i=1,\[Sigma]=1,orb=1 you get 0; if you set i=L+1,\[Sigma]=f,orb=Norb you get the total number of fermions in the state *)
CountFermions = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}
	},
	Total@Take[ (* sum *)
		Flatten[IntegerDigits[#,2,L]&/@state] (* flattened version of the binary representation of the state *)
		,L*(f*(orb-1)+\[Sigma]-1)+(i-1) (* sum up to this index in the flattened version of the state *)
	],
	CompilationTarget->"C"
];

(* sign accumulated by moving c_i_orb_\[Sigma] to the correct position when applying c_i_orb_\[Sigma]|state> or cdg_i_orb_\[Sigma]|state> *)
CSign = Compile[{
	{L,_Integer},{f,_Integer},{i,_Integer},{\[Sigma],_Integer},{orb,_Integer},{state,_Integer,1}
	},
	(-1)^Total@Take[ (* sum *)
		Flatten[IntegerDigits[#,2,L]&/@state] (* flattened version of the binary representation of the state *)
		,L*(f*(orb-1)+\[Sigma]-1)+(i-1) (* sum up to this index in the flattened version of the state *)
	],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True, RuntimeOptions->"Speed"
];

(* sign accumulated by moving c_i1_orb1_\[Sigma]1, c_i2_orb2_\[Sigma]2, c_i3_orb3_\[Sigma]3, ... to the correct positions when applying c_i1_orb1_\[Sigma]1 c_i2_orb2_\[Sigma]2 ...|state> or similar operators *)
CCSign = Compile[{
	{L,_Integer},{f,_Integer},{iList,_Integer,1},{\[Sigma]List,_Integer,1},{orbList,_Integer,1},{state,_Integer,1}
	},
	With[{
		binarystate = Flatten[IntegerDigits[#,2,L]&/@state],
		indexes = L*(f*(orbList-1)+\[Sigma]List-1)+(iList-1)
	},
	(-1)^Total[
		Total@Take[binarystate, #]&/@indexes
	]
	], CompilationTarget->"C", RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"
];

CCCCSign = CCSign;


(*               HOPPING FUNCTIONS             *)
(* gives True if hopping from (j, \[Sigma]2, orb2) to (i,\[Sigma]1,orb1) is possible, False otherwise *)
HopQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {j,_Integer}, {\[Sigma]1,_Integer}, {\[Sigma]2,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+\[Sigma]2]])[[j]]==1 && (IntegerDigits[#,2,L]&@state[[f*(orb1-1)+\[Sigma]1]])[[i]]==0, 
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True, RuntimeOptions->"Speed"
];

(* select states for which hopping (j, \[Sigma]2, orb2) \[Rule] (i, \[Sigma]1, orb1) is possible *)
HopSelect[L_, f_, i_, j_, \[Sigma]1_, \[Sigma]2_, orb1_, orb2_, stateList_] := Module[
	{criteria = HopQ[L,f,i,j,\[Sigma]1,\[Sigma]2,orb1,orb2,stateList]},
	Pick[stateList, criteria]
];

(* apply hopping operator to a given state *)
Hop = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {j,_Integer}, {\[Sigma]1,_Integer}, {\[Sigma]2,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}},
	MapAt[
		BitSet[#, L-i]&,
		MapAt[
			BitClear[#, L-j]&,
			state,
			f*(orb2-1)+\[Sigma]2
		],
		f*(orb1-1)+\[Sigma]1
	],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"
];

(* number of particles on site i with spin \[Sigma] in orbital orb *)
n = Compile[{
	{L,_Integer}, {f,_Integer}, {Norb,_Integer}, {i,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}},
	(IntegerDigits[#,2,L]&@state[[f*(orb-1)+\[Sigma]]])[[i]],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True, RuntimeOptions->"Speed"
];


(*           PAIR CREATION / ANNIHILATION FUNCTIONS          *)
(* gives True if it is possible to create a pair (i,orb1,\[Sigma]1) and (j,orb2,\[Sigma]2) and False otherwise *)
CreatePairQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {j,_Integer}, {\[Sigma]1,_Integer}, {\[Sigma]2,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+\[Sigma]1]])[[i]] == 0 && (IntegerDigits[#,2,L]&@state[[f*(orb2-1)+\[Sigma]2]])[[j]] == 0,
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True, RuntimeOptions->"Speed"
];

(* selects states for which pair creation (i,orb1,\[Sigma]1), (j,orb2,\[Sigma]2) is possible*)
CreatePairSelect[L_, f_, i_, j_, \[Sigma]1_, \[Sigma]2_, orb1_, orb2_, stateList_] := Module[
	{criteria = CreatePairQ[L, f, i, j, \[Sigma]1, \[Sigma]2, orb1, orb2, stateList]},
	Pick[stateList, criteria]
];

(* create a pair of particles (i,orb1,\[Sigma]1) (j,orb2,\[Sigma]2) and return the integer version of the states. *)
CreatePair = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {j,_Integer}, {\[Sigma]1,_Integer}, {\[Sigma]2,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	MapAt[
		BitSet[#, L-j]&,
		MapAt[
			BitSet[#, L-i]&,
			state,
			f*(orb1-1)+\[Sigma]1
		],
		f*(orb2-1)+\[Sigma]2
	],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"
];


(*          PAIR HOPPING FUNCTIONS         *)
(* gives True if it is possible to make a pair hopping from (i,orb2) to (i,orb1), False otherwise  *)
PairHoppingQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer},{orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+1]])[[i]] == 0 &&
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+2]])[[i]] == 0 &&
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+1]])[[i]] == 1 &&
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+2]])[[i]] == 1,
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True, RuntimeOptions->"Speed"
];

(* selects states for which pair hopping from (i,orb2) to (i,orb1) is possible*)
PairHoppingSelect[L_, f_, i_, orb1_, orb2_, stateList_] := Module[
	{criteria = PairHoppingQ[L,f,i,orb1,orb2,stateList]},
	Pick[stateList, criteria]
];

(* returns the integer form of the state obtained by pair hopping *)
PairHopping = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	MapAt[
		BitSet[#, L-i]&,
		MapAt[
			BitSet[#, L-i]&,
			MapAt[
				BitClear[#, L-i]&,
				MapAt[
					BitClear[#, L-i]&,
					state,
					f*(orb2-1)+2
				],
				f*(orb2-1)+1
			],
			f*(orb1-1)+1
		],
		f*(orb1-1)+2
	],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"
];


(*           SPIN EXCHANGE FUNCTIONS         *)
(* gives True if spin exchange on site i is possible, False otherwise: spin exchange = cdg_orb1_up c_orb1_dw cdg_orb2_dw c_orb2_up *)
SpinExchangeQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+1]])[[i]] == 1 && 
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+2]])[[i]] == 1 &&
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+2]])[[i]] == 0 && 
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+1]])[[i]] == 0, 
	CompilationTarget->"C", RuntimeAttributes->{Listable}, Parallelization->True, RuntimeOptions->"Speed"
];

(* select states for which spin exchange on site i is possible: spin exchange = cdg_orb1_up c_orb1_dw cdg_orb2_dw c_orb2_up *)
SpinExchangeSelect[L_, f_, i_, orb1_, orb2_, stateList_] := Module[
	{criteria = SpinExchangeQ[L,f,i,orb1,orb2,stateList]},
	Pick[stateList, criteria]
];

(* apply spin exchange operator to a given state *)
SpinExchange = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}},
	MapAt[
		BitSet[#, L-i]&,
		MapAt[
			BitSet[#, L-i]&,
			MapAt[
				BitClear[#, L-i]&,
				MapAt[
					BitClear[#, L-i]&,
					state,
					f*(orb2-1)+1
				],
				f*(orb1-1)+2
			],
			f*(orb1-1)+1
		],
		f*(orb2-1)+2
	],
	CompilationTarget->"C", RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"
];

(* gives right neighbor of a site j in a closed (open) chain of L elements *)
Neighbor = Compile[{
	{L,_Integer}, {j,_Integer}
	},
	Mod[j,L]+1,
	CompilationTarget->"C"
];


(*          BUILD THE IMPURITY HAMILTONIAN           *)
(* Non-local Hamiltonian blocks *)
HNonlocal[L_, f_, Norb_, Sectors_, EdMode_, OptionsPattern[]] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H = {};
	Do[
		Hsector = {};
		dim = Length[\[Psi]];
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Hblock = SparseArray[{}, {dim,dim}];
			Which[
				flag == "Bath",
				num = n[L, f, Norb, j, \[Sigma], orb, \[Psi]];(*local density*)
				Hblock += SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Hopping",
				\[Psi]1 = HopSelect[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]];
				If[Length[\[Psi]1] != 0, 
					\[Chi] = Hop[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]1];
					rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
					\[CapitalSigma] = CCSign[L, f, {1,j}, {\[Sigma],\[Sigma]}, {orb,orb}, \[Psi]1];
					Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					Hblock = Hblock + Hblock\[ConjugateTranspose];
				];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Superc" && (EdMode == "Superc" || EdMode == "FullSuperc"),
				If[\[Sigma] == 1,
					\[Psi]1 = CreatePairSelect[L, f, j, j, 1, 2, orb, orb, \[Psi]];
					If[Length[\[Psi]1] != 0, 
						\[Chi] = CreatePair[L, f, j, j, 1, 2, orb, orb, \[Psi]1];
						rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						\[CapitalSigma] = CCSign[L, f, {j,j}, {1,2}, {orb,orb}, \[Psi]1];
						Hblock += SparseArray[pos->\[CapitalSigma], {dim,dim}];
						Hblock = Hblock + Hblock\[ConjugateTranspose];
					];
					AppendTo[Hsector, Hblock];
				];,
			(* --------------------------------------------------------------- *)
				flag == "InterorbSuperc" && (EdMode == "InterorbSuperc" || EdMode == "FullSuperc"),
				If[\[Sigma] == 1,
				Do[
					If[
						orb2 > orb,
						\[Psi]1 = CreatePairSelect[L, f, j, j, 1, 2, orb, orb2, \[Psi]];
						If[Length[\[Psi]1] != 0,
							\[Chi] = CreatePair[L, f, j, j, 1, 2, orb, orb2, \[Psi]1];
							rows = \[Chi]/.dispatch;(* *)cols = \[Psi]1/.dispatch;(* *)pos = {rows,cols}\[Transpose];
							\[CapitalSigma] = CCSign[L, f, {j,j}, {1,2}, {orb,orb2}, \[Psi]1];
							Hblock += SparseArray[pos -> \[CapitalSigma], {dim,dim}];
						];
					]
				, {orb2,1,Norb}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];
			];
			];
		,{flag, {"Bath","Hopping","Superc","InterorbSuperc"}}, {orb,1,Norb}, {\[Sigma],1,f}, {j,OptionValue[Nimp]+1,L}];
		AppendTo[H, Hsector];
	,{\[Psi], Sectors}];
	H
];
Options[HNonlocal] = {Nimp -> 1};

(* prints useful info about the order of Hamiltonian blocks *)
HNonlocalInfo[L_, f_, Norb_, EdMode_] := Module[{},
	Print["The output blocks have the following order:"]
	Do[
		Print[flag,",  j=",j,",   \[Sigma]=",\[Sigma],",  orb=",orb]
	,{flag,{"Bath","Hopping"}}, {orb,1,Norb}, {\[Sigma],1,f}, {j,2,L}];
	If[
		EdMode == "Superc" || EdMode == "FullSuperc",
		Do[
			Print["Superc,  j=",j,",  orb=",orb]
		,{orb,1,Norb},{j,2,L}];
	];
	If[
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		Do[
			Do[
				If[
					orb2>orb,
					Print["InterorbSuperc,  j=",j,",  orbs=",orb," ",Mod[orb,Norb]+1]
				]
			,{orb2,1,Norb}]
		,{orb,1,Norb},{j,2,L}];
	];
];


(* Local Hamiltonian blocks *)
HLocal[L_, f_, Norb_, Sectors_, EdMode_, OptionsPattern[]] := Module[
	{H,Hsector,Hblock,rules,dispatch,num,dim,\[Psi]1,\[Chi],rows,cols,pos,\[CapitalSigma]},
	H = {};
	Do[
		Hsector = {};
		dim = Length[\[Psi]];
		Hsector = {};
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Which[
				flag == "Hubbard",
				Do[
					Hblock = SparseArray[{},{dim,dim}];
					Do[
						num = Sum[
							n[L, f, Norb, j, \[Sigma], orb, \[Psi]]
						, {\[Sigma], 1, f}];
						Hblock += SparseArray@DiagonalMatrix[0.5*num*(num-1)];
					, {j, 1, OptionValue[Nimp]}];
					AppendTo[Hsector, Hblock];
				, {orb, 1, Norb}],
			(* ---------------------------------- *)
				flag == "Interorb_Hubbard_Opposite_Spin" && Norb > 1,
				num = Sum[
					If[orbA != orbB,
						n[L,f,Norb,j,1,orbA,\[Psi]]*n[L,f,Norb,j,2,orbB,\[Psi]],
					(*else*)
						0]
				,{orbA,1,Norb}, {orbB,1,Norb}, {j,1,OptionValue[Nimp]}];
				Hblock = SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];,
			(* ---------------------------------- *)
				flag == "Interorb_Hubbard_Same_Spin" && Norb > 1,
				num = Sum[
					n[L,f,Norb,j,\[Sigma],1,\[Psi]]*n[L,f,Norb,j,\[Sigma],2,\[Psi]]
				,{\[Sigma],1,f}, {j,1,OptionValue[Nimp]}];
				Hblock = SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector,Hblock];,
			(* ---------------------------------- *)
				flag == "Pair_Hopping" && Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "Superc" || EdMode == "FullSuperc"),
				Hblock = SparseArray[{},{dim,dim}];
				Do[
					If[orb2 > orb1,
						\[Psi]1 = PairHoppingSelect[L, f, 1, orb1, orb2, \[Psi]];
						If[Length[\[Psi]1] == 0, Continue[];];
						\[Chi] = PairHopping[L, f, 1, orb1, orb2, \[Psi]1];
						rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						\[CapitalSigma] = CCSign[L, f, {j,j,j,j}, {1,2,1,2}, {orb1,orb1,orb2,orb2}, \[Psi]1];
						Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					];
				,{orb1,1,Norb}, {orb2,1,Norb}, {j,1,OptionValue[Nimp]}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector,Hblock];,
			(* ----------------------------------- *)
				flag == "Spin_Exchange" && Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc"),
				Hblock = SparseArray[{},{dim,dim}];
				Do[
					If[orb2 > orb1,
						\[Psi]1 = SpinExchangeSelect[L, f, 1, orb1, orb2, \[Psi]];
						If[Length[\[Psi]1] == 0, Continue[];];
						\[Chi] = SpinExchange[L, f, 1, orb1, orb2, \[Psi]1];
						rows = \[Chi]/.dispatch;(* *)cols = \[Psi]1/.dispatch;(* *)pos = {rows,cols}\[Transpose];
						\[CapitalSigma] = CCSign[L, f, {j,j,j,j}, {1,2,1,2}, {orb1,orb1,orb2,orb2}, \[Psi]1];
						Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					];
				,{orb1,1,Norb}, {orb2,1,Norb}, {j,1,OptionValue[Nimp]}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];,
			(* ---------------------------------- *)
				flag == "Energy_Shift",
				num = Sum[
					n[L,f,Norb,j,\[Sigma],orb,\[Psi]]
				,{\[Sigma],1,f}, {orb,1,Norb}, {j,1,OptionValue[Nimp]}];
				Hblock = SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];
			];
		,{flag, {"Hubbard","Interorb_Hubbard_Opposite_Spin","Interorb_Hubbard_Same_Spin","Pair_Hopping","Spin_Exchange","Energy_Shift"}}];
		AppendTo[H, Hsector];
	,{\[Psi],Sectors}];
	H
];
Options[HLocal] = {Nimp -> 1};


(* Get the Hamiltonian structure once for all *)
GetHamiltonian[L_, f_, Norb_, Sectors_, LoadHamiltonianQ_, HnonlocFile_, HlocFile_, EdMode_] := Module[
	{HnonlocBlocks, HlocBlocks},
	If[LoadHamiltonianQ,
		Print["Getting Hamiltonians from file"];
		HnonlocBlocks = Import[HnonlocFile];
		HlocBlocks = Import[HlocFile];
		Print["Done! Let's get started! \n"],
	(* else *)
		Print["Computing Hamiltonians..."];
		Print["Time: ", First @ AbsoluteTiming[
			HnonlocBlocks = HNonlocal[L, f, Norb, Sectors, EdMode];
			Export[HnonlocFile, HnonlocBlocks];
			HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode];
			Export[HlocFile, HlocBlocks];
		]," sec."];
		Print["Done! Let's get started! \n"];
	];
	{HnonlocBlocks, HlocBlocks}
];

(* build the impurity Hamiltonian from the local and nonlocal blocks and the respective parameters *)
HImp[Norb_, HnonlocBlocks_, HlocBlocks_, BathParameters_, InteractionParameters_, EdMode_] := Module[
	{EffectiveInteractionParameters, FlatBathParameters = Flatten[BathParameters], Hloc, Hnonloc},
	EffectiveInteractionParameters = Which[
		Norb == 1,
		(* with 1 orbital, delete Jse, Jph, Usec, Ust, at positions -2, -3, -4, -5 *)
		Delete[InteractionParameters, {{-2},{-3},{-4},{-5}}],
	(* --------------------------------------------------- *)
		Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc"),
		(* this is the most complicated scenario, you don't remove anything here *)
		InteractionParameters,
	(* --------------------------------------------------- *)
		Norb > 1 && EdMode == "Superc", 
		(* pair hopping is possible, but spin exchange is not *)
		Delete[InteractionParameters, -2],
	(* --------------------------------------------------- *)
		Norb > 1 && EdMode == "Normal",
		(* pair hopping and spin exchange are not possible *)
		Delete[InteractionParameters, {{-2},{-3}}]
	];
	(* build non local Hamiltonian -> avoid Dot[] as it returns dense array *)
	Hnonloc = SparseArray[#]&/@(
		Sum[
			FlatBathParameters[[i]]*#[[i]],
		{i, 1, Length@FlatBathParameters}]&/@HnonlocBlocks
	);
	(* build local Hamiltonian -> avoid Dot[] as it returns dense array *)
	Hloc = SparseArray[#]&/@(
		Sum[
			EffectiveInteractionParameters[[i]]*#[[i]],
		{i, 1, Length@EffectiveInteractionParameters}]&/@HlocBlocks
	);
	Hnonloc + Hloc
];


(*       BUILD CHAIN HAMILTONIAN        *)
Hnonint[L_, f_, Norb_, Sectors_, EdMode_, OptionsPattern[]] := Module[
	{\[Psi]1,\[Chi],H,Hblock,Hsector,dim,rules,dispatch,cols,rows,pos,\[CapitalSigma],num},
	H = {};
	Do[
		Hsector = {};
		dim = Length[\[Psi]];
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Hblock = SparseArray[{}, {dim,dim}];
			Which[
				flag == "Potential",
				num = n[L, f, Norb, j, \[Sigma], orb, \[Psi]];(*local density*)
				Hblock += SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Hopping",
				If[j == L && !OptionValue[RealPBC], Continue[];];
				\[Psi]1 = HopSelect[L, f, j, Neighbor[L,j], \[Sigma], \[Sigma], orb, orb, \[Psi]];
				If[Length[\[Psi]1] == 0, AppendTo[Hsector, Hblock]; Continue[];];
				\[Chi] = Hop[L, f, j, Neighbor[L,j], \[Sigma], \[Sigma], orb, orb, \[Psi]1];
				rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
				\[CapitalSigma] = CCSign[L,f,{j,Neighbor[L,j]},{\[Sigma],\[Sigma]},{orb,orb},\[Psi]1];
				If[j == L, \[CapitalSigma] = -\[CapitalSigma]];(* if j=L, we are applying cdg_L c_1. When moving cdg_L before position L, it jumps over c_1 and the sign changes! *)
				If[
					OptionValue[RealPhase] == {0},
					Hblock += SparseArray[pos -> \[CapitalSigma], {dim,dim}];,
				(* else *)
					Hblock += SparseArray[pos -> \[CapitalSigma]*Exp[I*OptionValue[RealPhase][[\[Sigma]]]], {dim,dim}]
				];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Raman" && EdMode == "Raman",
				If[\[Sigma] == f && !OptionValue[SyntheticPBC], Continue[];];
				\[Psi]1 = HopSelect[L, f, j, j, \[Sigma], Neighbor[f,\[Sigma]], orb, orb, \[Psi]];
				If[Length[\[Psi]1] == 0, AppendTo[Hsector, Hblock]; Continue[];];
				\[Chi] = Hop[L, f, j, j, \[Sigma], Neighbor[f,\[Sigma]], orb, orb, \[Psi]1];
				rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
				\[CapitalSigma] = CCSign[L,f,{j,j},{\[Sigma],Neighbor[f,\[Sigma]]},{orb,orb},\[Psi]1];
				If[\[Sigma] == f, \[CapitalSigma] = -\[CapitalSigma]];(* if \[Sigma]=f, we are applying cdg_f c_1. When moving cdg_f before position f, it jumps over c_1 and the sign changes! *)
				Hblock += SparseArray[pos -> \[CapitalSigma]*Exp[I*OptionValue[SyntheticPhase]*j], {dim,dim}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];	
			];
		,{flag, {"Potential","Hopping","Raman"}}, {orb,1,Norb}, {\[Sigma],1,f}, {j,1,L}];
		AppendTo[H, Hsector];
	,{\[Psi], Sectors}];
	H
];
Options[Hnonint] = {RealPBC -> True, SyntheticPBC -> False, RealPhase -> {0}, SyntheticPhase -> 0};


(*                       LANCZOS                       *)
Lanczos[H_, StartingVector_, OptionsPattern[]] := Module[{
	miniter = OptionValue[MinIter],
	maxiter = OptionValue[MaxIter],
	\[Epsilon] = OptionValue[ConvergenceThreshold],
	shift = OptionValue[Shift],
	a,b,a0,b1,v,w,HKrilov,E0old,E0new,nfinal
	},
	(* initialize array of a_n :  a[1]=a_0 , a[1+n]=a_n *)
	a = ConstantArray[0, maxiter+1];
	(*initialize array of b_n : b[n]=b_n *)
	b = ConstantArray[0, maxiter];
	(* initialize starting vector *)
	v = Normalize[StartingVector];
	(* initialize a new vector w = Hv *)
	w = H . v;
	a0 = (Conjugate[v]) . w;
	w = w - a0*v;
	b1 = Norm[w];
	a[[1]] = a0; b[[1]] = b1;
	E0old = a0;
	nfinal = maxiter;
	Do[
		If[b[[n]] < \[Epsilon] && n >= miniter,
			nfinal = n; Break[];
		];
		w = w/b[[n]];(* w=Subscript[v, n] *)
		v = -b[[n]]*v;(* v=-Subscript[b, n]Subscript[v, n-1]*)
		{v,w} = {w,v};
		w = w + H . v;(* w=Subscript[Hv, n]-Subscript[b, n]Subscript[v, n-1] *)
		a[[n+1]] = (Conjugate@v) . w; (*Subscript[a, n] = Subscript[v, n]Subscript[Hv, n]-Subscript[b, n]Subscript[v, n]Subscript[v, n-1]*)
		w = w - a[[n+1]] * v;
		If[n < maxiter,
			b[[n+1]] = Norm[w];
		];
		(*Build hamiltonian in the Krylov subspace*)
		HKrilov = SparseArray[{
			Band[{1,2}] -> Take[b, n],
			Band[{2,1}] -> Take[b, n],
			Band[{1,1}] -> Take[a, n+1]
		},{n+1, n+1}];
		E0new = If[shift==0,
			Min@Eigenvalues@HKrilov,
		(*else*)	
			Eigenvalues[HKrilov-DiagonalMatrix[ConstantArray[shift,1+n]]][[1]]+shift
		];
		If[Abs[E0new-E0old] < \[Epsilon] && n > miniter,
			nfinal = n; Break[];
		];
		E0old=E0new;
	,{n,1,maxiter}];
	a = Take[a, nfinal+1];
	b = Take[b, nfinal];
	{E0new,a,b}
];
Options[Lanczos] = {ConvergenceThreshold -> 1.0*10^(-8), MinIter -> 2, MaxIter -> 2000, Shift -> 0};

(* i,j element of the inverse matrix *)
InverseElement[m_, {i_,j_}] := (-1)^(i+j)Det[Drop[m,{j},{i}]]/Det[m];


(* ANALYTIC EVALUATION OF NONINTERACTING IMPURITY GREEN FUNCTION *)
(* non-interacting Green function *)
GreenFunction0[L_, f_,\[Mu]_, symbols_, z_, EdMode_] := Module[
	{e, V, \[CapitalDelta], H, G},
	Which[
		EdMode == "Normal",
		e = symbols[[1;;L-1]];
		V = symbols[[L;;2(L-1)]];
		(* the spinor is (d, c_1, c_2, ...) *)
		H = SparseArray[{
			{i_,j_}/;(i==1&&j>1):>V[[j-1]]
		},
		{L,L}];
		H = H + H\[Transpose];
		H += SparseArray[{
			{1,1}->-\[Mu],
		{i_,i_}/;(i>1):>e[[i-1]]
		},
		{L,L}];
		G = InverseElement[SparseArray[z*IdentityMatrix[L]-H], {1,1}];,
	(* ------------------------ *)
		EdMode == "Superc",
		e = symbols[[1;;L-1]];
		V = symbols[[L;;2(L-1)]];
		\[CapitalDelta] = symbols[[2L-1;;3(L-1)]];
		(* the spinor is (d_up, ddg_dw, c_1up, cdg_1dw, c_2up, cdg_2dw, ... *)
		H = SparseArray[{
			{i_,j_}/;(j==i+1&&i>2&&Mod[i,f]==1):>\[CapitalDelta][[Quotient[i,f]]],
			{i_,j_}/;(i==1&&j>2&&Mod[j,f]==1):>V[[Quotient[j,f]]],
			{i_,j_}/;(i==2&&j>2&&Mod[j,f]==0):>-V[[Quotient[j-1,f]]]
		},
		{f*L, f*L}];
		H = H + H\[Transpose];
		H += SparseArray[{
			{1, 1}->-\[Mu], {2,2}->\[Mu],
			{i_,i_}/;(i>2&&Mod[i,f]==1):>e[[Quotient[i-1,f]]],
			{i_,i_}/;(i>2&&Mod[i,f]==0):>-e[[Quotient[(i-2),f]]]
		},
		{f*L,f*L}
		];
		G = Table[
			InverseElement[SparseArray[z*IdentityMatrix[f*L] - H], {i,j}]
		,{i,1,2},{j,1,2}];(* the IMPURITY part of the Green function is the 2x2 top left block *)
	];
G
];

(* Weiss field: non-interacting inverse of the Green function *)
WeissField[L_, f_, \[Mu]_, symbols_, z_, EdMode_] := With[
	{G0 = GreenFunction0[L, f, \[Mu], symbols, z, EdMode]},
	Which[
		EdMode == "Normal", FullSimplify[1/G0],
		EdMode == "Superc", FullSimplify[Inverse[G0]]
	]
];


(*                APPLY CDG / C TO STATES             *)
(* find the quantum number of the final state after application of operator_{j,\[Sigma],orb} to a state living in the sector qns *)
FinalSector[L_, f_, Norb_, j_, \[Sigma]_, orb_, qns_, operator_, EdMode_] := Module[
	{newqns = qns},
	Which[
		operator == "Creation",
		Which[
			EdMode == "Normal",
			If[qns[[f*(orb-1)+\[Sigma]]] == L, Return[]];  (* trivial case: return Null if final sector does not exist *)
			newqns[[f*(orb-1)+\[Sigma]]] += 1;,
		(* ---------------------------- *)
			EdMode == "InterorbNormal",
			If[qns[[\[Sigma]]] == Norb*L, Return[]]; (* trivial case: return Null if sector does not exist *)
			newqns[[\[Sigma]]] += 1;,
		(* ---------------------------- *)
			EdMode == "Superc",
			If[f>2, Return["error. f>2 not supported with EdMode = ''Superc''"];];
			If[(qns[[orb]] == -L && \[Sigma] == 2) || (qns[[orb]] == L && \[Sigma] == 1), Return[]];  (* trivial case: return Null if sector does not exist *)
			Which[
				\[Sigma]==1, newqns[[orb]] += 1,
				\[Sigma]==2, newqns[[orb]] -= 1
			];,
		(* ---------------------------- *)
			EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
			If[f>2, Return["error. f>2 not supported with EdMode = ''FullSuperc'' or ''InterorbSuperc''"];];
			If[(qns == -Norb*L && \[Sigma]==2) || (qns == Norb*L && \[Sigma]==1), Return[]];  (* trivial case: return Null if final sector does not exist. *)
			Which[
				\[Sigma]==1, newqns+=1,
				\[Sigma]==2, newqns-=1
			];
		],
	(* ------------------------------------------ *)	
		operator == "Annihilation",
		Which[
			EdMode == "Normal",
			If[qns[[f*(orb-1)+\[Sigma]]] == 0, Return[]];  (* trivial case *)
			newqns[[f*(orb-1)+\[Sigma]]] -= 1;,
		(* ---------------------------- *)
			EdMode == "InterorbNormal",
			If[qns[[\[Sigma]]] == 0, Return[]];  (* trivial case *)
			newqns[[\[Sigma]]] -= 1;,
		(* ---------------------------- *)
			EdMode == "Superc",
			If[f>2, Return["error. f>2 not supported with EdMode = ''Superc''"];];
			If[(qns[[orb]] == -L && \[Sigma] == 1) || (qns[[orb]] == L && \[Sigma] == 2), Return[]];  (* trivial case *)
			Which[
				\[Sigma]==1, newqns[[orb]] -= 1,
				\[Sigma]==2, newqns[[orb]] += 1
			];,
		(* ---------------------------- *)
			EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
			If[f>2, Return["error. f>2 not supported with EdMode = ''FullSuperc'' or ''InterorbSuperc''"];];
			If[(qns == -Norb*L && \[Sigma]==1) || (qns == Norb*L && \[Sigma]==2), Return[]];  (* trivial case *)
			Which[
				\[Sigma]==1, newqns -= 1,
				\[Sigma]==2, newqns += 1
			];
		],
	(* ------------------------------------------ *)	
		operator == "Density",
		newqns = qns;	
	];
	newqns
];

(* apply cdg_{j,\[Sigma],orb} |gs>, where |gs> belongs to the sector with quantum numbers qns and give the resulting vector resized to fit the dimension of the sector obtained adding a particle with state label (j,\[Sigma],orb)*)
ApplyCdg[L_, f_, Norb_, j_, \[Sigma]_, orb_, gs_, qns_, Sectors_, SectorsDispatch_, EdMode_] := Module[
	{newqns = qns,startingsector = Sectors[[qns/.SectorsDispatch]],finalsector, newdim,sign,dispatch,pos,newpos,coeff,\[Psi]1,\[Chi]},
	(* check which states of the starting sector can host the extra particle *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,startingsector],1]];
	\[Psi]1 = CreateParticleSelect[L, f, j, \[Sigma], orb, startingsector];
	pos = \[Psi]1/.dispatch;
	(* compute the correct signs obtained moving cdg to the correct position  *)
	sign = CSign[L, f, j, \[Sigma], orb, \[Psi]1];
	(* list of coefficients that remain non vanishing *)
	coeff = gs[[pos]]*sign;
	(* build the final sector *)
	Which[
		EdMode=="Normal",
		If[qns[[f*(orb-1)+\[Sigma]]]==L,Return[0]];  (* trivial case *)
		newqns[[f*(orb-1)+\[Sigma]]]+=1;,
	(* ---------------------------- *)
		EdMode == "InterorbNormal",
		If[qns[[\[Sigma]]] == Norb*L, Return[0]];  (* trivial case *)
		newqns[[\[Sigma]]]+=1;,
	(* ---------------------------- *)
		EdMode == "Superc",
		If[f>2, Return["error. f>2 not supported with EdMode = ''Superc''"];];
		If[(qns[[orb]] == -L && \[Sigma] == 2) || (qns[[orb]] == L && \[Sigma] == 1),Return[0]];  (* trivial case *)
		Which[
			\[Sigma]==1, newqns[[orb]]+=1,
			\[Sigma]==2, newqns[[orb]]-=1
		],
	(* ---------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		If[f>2, Return["error. f>2 not supported with EdMode = ''FullSuperc'' or ''InterorbSuperc''"];];
		If[(qns==-Norb*L&&\[Sigma]==2)||(qns==Norb*L&&\[Sigma]==1),Return[0]];  (* trivial case *)
		Which[
			\[Sigma]==1, newqns+=1,
			\[Sigma]==2, newqns-=1
		]
	];
	finalsector = Sectors[[newqns/.SectorsDispatch]];
	newdim = Length[finalsector];
	(* create a dispatch that labels all these states *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1]];(* find the new positions where the entries of gs should go after applying cdg_spin *)
	\[Chi] = cdg[L,f,j,\[Sigma],orb,\[Psi]1];
	newpos = \[Chi]/.dispatch;
	(* create a list of rules and define the resulting array *)
	SparseArray[Thread[newpos->coeff],newdim]
];

(* apply c_{j,\[Sigma],orb} |gs>, where |gs> belongs to the sector with quantum numbers qns and give the resulting vector resized to fit the dimension of the sector obtained removing a particle with state label (j,\[Sigma],orb)*)
ApplyC[L_, f_, Norb_, j_, \[Sigma]_, orb_, gs_, qns_, Sectors_, SectorsDispatch_, EdMode_] := Module[
	{newqns = qns,startingsector = Sectors[[qns/.SectorsDispatch]],finalsector, newdim,sign,dispatch,pos,newpos,coeff,\[Psi]1,\[Chi]},
	(* check which states of the starting sector can host the extra particle *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,startingsector],1]];
	\[Psi]1 = DestroyParticleSelect[L, f, j, \[Sigma], orb, startingsector];
	pos = \[Psi]1/.dispatch;
	(* compute the correct signs obtained moving cdg to the correct position  *)
	sign = CSign[L,f,j,\[Sigma],orb,\[Psi]1];
	(* list of coefficients that remain non vanishing *)
	coeff = gs[[pos]]*sign;
	(* build the final sector *)
	Which[
		EdMode=="Normal",
		If[qns[[f*(orb-1)+\[Sigma]]]==0, Return[0]];  (* trivial case *)
		newqns[[f*(orb-1)+\[Sigma]]]-=1;,
	(* ---------------------------- *)
		EdMode =="InterorbNormal",
		If[qns[[\[Sigma]]]==0, Return[0]];  (* trivial case *)
		newqns[[\[Sigma]]]-=1;,
	(* ---------------------------- *)
		EdMode == "Superc",
		If[f>2, Return["error. f>2 not supported with EdMode = ''Superc''"];];
		If[(qns[[orb]]==-L&&\[Sigma]==1)||(qns[[orb]]==L&&\[Sigma]==2), Return[0]];  (* trivial case *)
		Which[
			\[Sigma]==1, newqns[[orb]]-=1,
			\[Sigma]==2, newqns[[orb]]+=1
		],
	(* ---------------------------- *)
		EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		If[f>2, Return["error. f>2 not supported with EdMode = ''FullSuperc'' or ''InterorbSuperc''"];];
		If[(qns==-Norb*L&&\[Sigma]==1)||(qns==Norb*L&&\[Sigma]==2),Return[0]];  (* trivial case *)
		Which[
			\[Sigma]==1, newqns-=1,
			\[Sigma]==2, newqns+=1
		]
	];
	finalsector = Sectors[[newqns/.SectorsDispatch]];
	newdim = Length[finalsector];
	(* create a dispatch that labels all these states *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1]];(* find the new positions where the entries of gs should go after applying cdg_spin *)
	\[Chi] = c[L,f,j,\[Sigma],orb,\[Psi]1];
	newpos = \[Chi]/.dispatch;
	(* create a list of rules and define the resulting array *)
	SparseArray[Thread[newpos->coeff], newdim]
];



(*      EXACT CALCULATION OF GREEN FUNCTION       *)
GreenFunctionED[L_, f_, Norb_, {i_,j_}, \[Sigma]_, orb_, Sectors_, QnsSectorList_, eigs_, T_, zlist_, EdMode_, OptionsPattern[]] := Module[
	{newqns, SectorsDispatch, startingsector, finalsector, dim, newdim, dispatch, sign, rows, cols, pos, \[Psi]1, \[Psi]2, \[Chi]1, \[Chi]2, A, B, P, Q, S, Z=1., e0, energies, eigenstates, G},
	(* get energies and eigenstates *)
	energies = eigs[[All, 1]];
	eigenstates = eigs[[All, 2]];
	e0 = Min[Flatten[energies]];
	(* get normalization factor if required *)
	If[OptionValue[NormalizedFunction], Z = Total[Exp[-(Flatten[energies]-e0)/T]]; ];
	(* get dispatch of the quantum numbers and initialize G.F. *)
	SectorsDispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,QnsSectorList],1]];
	G = ConstantArray[0, Length[zlist]];
	(* loop over all the sectors *)
	Do[
		(* get starting sector *)
		startingsector = Sectors[[qns/.SectorsDispatch]];
		dim = Length[startingsector];
		dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&, startingsector],1]];
		(* select those states for which you can create a particle with index j,\[Sigma],orb *)
		\[Psi]1 = CreateParticleSelect[L, f, j, \[Sigma], orb, startingsector];
		cols = \[Psi]1/.dispatch;
		sign = CSign[L, f, j, \[Sigma], orb, \[Psi]1];
		(* get final sector *)
		newqns = FinalSector[L, f, Norb, 1, \[Sigma], orb, qns, "Creation", EdMode];
		If[newqns == Null, Continue[];];
		finalsector = Sectors[[newqns/.SectorsDispatch]];
		newdim = Length[finalsector];
		dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1]];
		\[Chi]1 = cdg[L,f,j,\[Sigma],orb,\[Psi]1];
		rows = \[Chi]1/.dispatch;
		(* build the matrix Smn = <n|c_i,\[Sigma]1,orb1|m><m|cdg_j,\[Sigma]2,orb2|n> * (e^-\[Beta]Em + e^-\[Beta]En)/(z + En - Em) *)
		(* where |n> is the n-th eigenstate that belongs to the sector with quantum number qns, and |m> is in a sector with one more particle. *)
		(* let Amn = <m|cdg_j,\[Sigma]2,orb2|n> and Bnm = <n|c_i,\[Sigma]1,orb1|m> *)
		pos = {rows, cols}\[Transpose];
		A = SparseArray[Thread[pos->sign], {newdim, dim}];(* A in the canonical basis *)
		P = Transpose[eigenstates[[qns/.SectorsDispatch]]];(* change of basis in starting sector *)
		Q = Transpose[eigenstates[[newqns/.SectorsDispatch]]];(* change of basis in final sector *)
		A = Q\[HermitianConjugate] . A . P;(* A in the eigenstates basis *)
	(* ------------------------------------- *)
		\[Psi]2 = DestroyParticleSelect[L, f, i, \[Sigma], orb, finalsector];
		rows = \[Psi]2/.dispatch;
		sign = CSign[L,f,i,\[Sigma],orb,\[Psi]2];
		dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,startingsector],1]];(* find the new positions where the entries of gs should go after applying cdg_spin *)
		\[Chi]2 = c[L,f,i,\[Sigma],orb,\[Psi]2];
		cols = \[Chi]2/.dispatch;
		pos = {rows, cols}\[Transpose];
		B = SparseArray[Thread[pos->sign], {newdim, dim}];(* B in the canonical basis *)
		B = Q\[HermitianConjugate] . B . P;
	(* ------------------------------------- *)
		S = A * B * Table[
			Exp[-(energies[[qns/.SectorsDispatch]][[n]]-e0)/T] + Exp[-(energies[[newqns/.SectorsDispatch]][[m]]-e0)/T]
		,{m,newdim},{n,dim}];
		G = G + (Total[#,2]&/@(
			S * Table[
				1./(# + energies[[qns/.SectorsDispatch]][[n]] - energies[[newqns/.SectorsDispatch]][[m]])
			,{m, newdim}, {n, dim}]&/@zlist));
	,{qns, QnsSectorList}];
	G/Z
];
Options[GreenFunctionED] = {NormalizedFunction -> False};

(* impurity Green function (spin and orbital diagonal part) *)
GreenFunctionImpurity[L_, f_, Norb_, \[Sigma]_, orb_, Egs_, gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_]:=Module[
	{startingsector = Sectors[[GsQns/.SectorsDispatch]],
	newqns = GsQns,
	GF = ConstantArray[0, Length[zlist]],
	finalsector, newdim, sign, dispatch, pos, newpos, coeff, \[Psi]1, \[Chi], cdggs, cgs, norm, H, E0, a, b},
(* STEP 1:   apply Cdg_\[Sigma],orb | gs > *)	
	(* 1.1:   check which states of the starting sector can host an extra particle *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,startingsector],1]];
	\[Psi]1 = CreateParticleSelect[L, f, 1, \[Sigma], orb, startingsector];
	pos = \[Psi]1/.dispatch;
	(* 1.2:   compute the resulting vector components *)
	sign = CSign[L, f, 1, \[Sigma], orb, \[Psi]1]; (* correct signs obtained moving cdg to the correct position *)
	coeff = gs[[pos]]*sign;(* list of coefficients that remain non vanishing *)
	(* 1.3:   build the final sector with an extra particle *)
	newqns = FinalSector[L, f, Norb, 1, \[Sigma], orb, GsQns, "Creation", EdMode];
	If[newqns == Null, Return[GF]];
	finalsector = Sectors[[newqns/.SectorsDispatch]];
	newdim = Length[finalsector];
	(* 1.4:  find the new positions where the entries of gs should go after applying cdg_spin *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1]];
	\[Chi] = cdg[L, f, 1, \[Sigma], orb, \[Psi]1];
	newpos = \[Chi]/.dispatch;
	(* 1.5:   resulting vector *)
	cdggs = SparseArray[Thread[newpos->coeff], newdim];
	norm = (Conjugate@cdggs) . cdggs;
(* STEP 2:    *)
	(* 2.1:   find Hamiltonian block in the final sector and the corresponding Krylov matrix *)
	H = Hsectors[[newqns/.SectorsDispatch]];
	{E0, a, b} = Lanczos[H, cdggs/norm];
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GF += norm*(
		InverseElement[
			SparseArray[(# + Egs) * IdentityMatrix[Length[a]] - H]
		, {1, 1}]&/@zlist);(* <gs|c cdg |gs> ((z - (HKrylov - Egs))^-1)[[1,1]] *)
(* STEP 3:   apply C_\[Sigma],orb | gs > *)	
	(* 3.1:   check which states of the starting sector can host an extra particle *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,startingsector],1]];
	\[Psi]1 = DestroyParticleSelect[L, f, 1, \[Sigma], orb, startingsector];
	pos = \[Psi]1/.dispatch;
	(* 3.2:   compute the resulting vector components *)
	sign = CSign[L, f, 1, \[Sigma], orb, \[Psi]1]; (* correct signs obtained moving c to the correct position *)
	coeff = gs[[pos]]*sign;(* list of coefficients that remain non vanishing *)
	(* 3.3:   build the final sector with an extra hole *)
	newqns = FinalSector[L, f, Norb, 1, \[Sigma], orb, GsQns, "Annihilation", EdMode];
	If[newqns == Null, Return[GF]]; 
	finalsector = Sectors[[newqns/.SectorsDispatch]];
	newdim = Length[finalsector];
	(* 3.4:  find the new positions where the entries of gs should go after applying cdg_spin *)
	dispatch = Dispatch[Flatten[MapIndexed[{#1->#2[[1]]}&,finalsector],1]];
	\[Chi] = c[L, f, 1, \[Sigma], orb, \[Psi]1];
	newpos = \[Chi]/.dispatch;
	(* 3.5:   resulting vector *)
	cgs = SparseArray[Thread[newpos->coeff], newdim];
	norm = (Conjugate@cgs) . cgs;
(* STEP 4:    *)
	(* 4.1:   find Hamiltonian block in the final sector and the corresponding Krylov matrix *)
	H = Hsectors[[newqns/.SectorsDispatch]];
	{E0, a, b} = Lanczos[H, cgs/norm];
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GF += norm*(
		InverseElement[
			SparseArray[(# - Egs) * IdentityMatrix[Length[a]] + H]
		, {1, 1}]&/@zlist);(* <gs|cdg c |gs> ((z + (HKrylov - Egs))^-1)[[1,1]] *)
	GF
];

(* compute the impurity Green function in the Nambu formalism *)
GreenFunctionImpurityNambu[L_, f_, Norb_, orb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_] := Module[{
	(* compute cdg_s|gs> and c_s|gs> for s = up, dw *)
	cdgup = ApplyCdg[L, f, Norb, 1, 1, orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode],
	cdgdw = ApplyCdg[L, f, Norb, 1, 2, orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode],
	cup = ApplyC[L, f, Norb, 1, 1, orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode],
	cdw = ApplyC[L, f, Norb, 1, 2, orb, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode],
	Odggs, Ogs, Pdggs, Pgs,
	newqns, H, E0, a, b,
	GFAparticle, GFAhole, GFBparticle, GFBhole, GFOparticle, GFOhole, GFO, GFPparticle, GFPhole, GFP, GFA, GFB, GF12, GF21},	
(* compute all the main contributions to the Green function *)
(*          G_O(z) "Particle" contribution             *)
	Odggs = cdgup + cdw;(* apply Odg|gs> = (Cdg_up + C_dw)|gs> *)
	newqns = FinalSector[L, f, Norb, 1, 1, orb, GsQns, "Creation", EdMode];
	H = Hsectors[[newqns/.SectorsDispatch]];(*Hamiltonian on that sector*)
	{E0,a,b} = Lanczos[H, Normalize[Odggs]];(* Apply Lanczos starting from Odg|gs> *)
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GFOparticle = (Norm[Odggs]^2)*(
		InverseElement[
			SparseArray[(# + Egs) * IdentityMatrix[Length[a]] - H]
		, {1, 1}] &/@ zlist);
(*           G_O(z) "Hole" contribution               *)
	Ogs = cup + cdgdw;(* apply (C_up + Cdg_dw)|gs> = O|gs> *)
	newqns = FinalSector[L, f, Norb, 1, 1, orb, GsQns, "Annihilation", EdMode];
	H = Hsectors[[newqns/.SectorsDispatch]];
	{E0,a,b} = Lanczos[H, Normalize[Ogs]];
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GFOhole = (Norm[Ogs]^2)*(
		InverseElement[
			SparseArray[(# - Egs) * IdentityMatrix[Length[a]] + H]
		, {1, 1}] &/@ zlist);
(*         G_P(z) "Particle" contribution           *)
	Pdggs = cdgup + I*cdw;(* apply (Cdg_up + I*C_dw)|gs> = Pdg|gs> *)
	newqns = FinalSector[L, f, Norb, 1, 1, orb, GsQns, "Creation", EdMode];(* if you create an up fermion or destroy a down fermion, you go from sz to sz+1 *)
	H = Hsectors[[newqns/.SectorsDispatch]];(* Hamiltonian on that sector *)
	{E0,a,b} = Lanczos[H, Normalize[Pdggs]];(*Apply Lanczos starting from Adg|gs> *)
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GFPparticle = (Norm[Pdggs]^2)*(
		InverseElement[
			SparseArray[(# + Egs) * IdentityMatrix[Length[a]] - H]
		, {1, 1}]&/@zlist);
(*          G_P(z) "Hole" contribution             *)
	Pgs = cup - I*cdgdw;(*apply (C_up - I*Cdg_dw)|gs> = P|gs> *)
	newqns = FinalSector[L, f, Norb, 1, 1, orb, GsQns, "Annihilation", EdMode];(*if you remove an up fermion or create a down fermion, you go from sz to sz-1*)
	H = Hsectors[[newqns/.SectorsDispatch]];
	{E0,a,b} = Lanczos[H, Normalize[Pgs]];
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GFPhole = (Norm[Pgs]^2)*(
		InverseElement[
			SparseArray[(# - Egs) * IdentityMatrix[Length[a]] + H]
		, {1, 1}]&/@zlist);
(*          G_up,up(z) "Particle" contribution             *)
	newqns = FinalSector[L, f, Norb, 1, 1, orb, GsQns, "Creation", EdMode];
	H = Hsectors[[newqns/.SectorsDispatch]];
	{E0,a,b} = Lanczos[H, Normalize[cdgup]];
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GFAparticle = (Norm[cdgup]^2)*(
		InverseElement[
			SparseArray[(# + Egs) * IdentityMatrix[Length[a]] - H]
		, {1, 1}]&/@zlist);
(*          G_up,up(z) "Hole" contribution             *)
	newqns = FinalSector[L, f, Norb, 1, 1, orb, GsQns, "Annihilation", EdMode];
	H = Hsectors[[newqns/.SectorsDispatch]];
	{E0,a,b} = Lanczos[H, Normalize[cup]];
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GFAhole = (Norm[cup]^2)*(
		InverseElement[
			SparseArray[(# - Egs) * IdentityMatrix[Length[a]] + H]
		, {1, 1}]&/@zlist);
(*          G_dw,dw(-z) "Particle" contribution             *)
	newqns = FinalSector[L, f, Norb, 1, 2, orb, GsQns, "Creation", EdMode];
	H = Hsectors[[newqns/.SectorsDispatch]];
	{E0,a,b} = Lanczos[H, Normalize[cdgdw]];
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GFBparticle = (Norm[cdgdw]^2)*(
		InverseElement[
			SparseArray[(# + Egs) * IdentityMatrix[Length[a]] - H]
		, {1, 1}]&/@(-zlist));		
(*          G_dw,dw(-z) "Hole" contribution             *)
	newqns = FinalSector[L, f, Norb, 1, 2, orb, GsQns, "Annihilation", EdMode];	
	H = Hsectors[[newqns/.SectorsDispatch]];
	{E0,a,b} = Lanczos[H, Normalize[cdw]];
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a]];(* Krylov matrix in the final sector *)
	GFBhole = (Norm[cdw]^2)*(
		InverseElement[
			SparseArray[(# - Egs) * IdentityMatrix[Length[a]] + H]
		, {1, 1}]&/@(-zlist));					
(*          Build Green function components               *)	
	GFA = GFAparticle + GFAhole;(* G_upup(z) *)
	GFB = -(GFBparticle + GFBhole);(* -G_dwdw(-z) *)
	GFO = GFOparticle + GFOhole;(* G_0(z) *)
	GFP = GFPparticle + GFPhole;(* G_P(z) *)
	(*compute the off-diagonal part using the diagonal part and GFO, GFP*)
	GF12 = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GFA + GFB));
	GF21 = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GFA + GFB));
	(* return a list with NMatsubara 2x2 matrices, each of them being the G.F. at that specific frequency *)
	Partition[#,2]&/@({GFA,GF12,GF21,GFB}\[Transpose])
];

InverseGreenFunction[L_, f_, Norb_, \[Sigma]_, orb_, Egs_, gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_] := With[{
	G = Which[
		EdMode == "Normal", GreenFunctionImpurity[L, f, Norb, \[Sigma], orb, Egs, gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist],
		EdMode == "Superc", GreenFunctionImpurityNambu[L, f, Norb, orb, Egs, gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist]
	]},
	Which[
		EdMode == "Normal", 1./G,
		EdMode == "Superc", Inverse/@G
	]
];

(* density of states of the infinite dimensional Bethe lattice *)
DoSBethe = Compile[
	{{\[Epsilon], _Real}, {DBethe, _Real}},
	(2./(Pi*DBethe^2))*Sqrt[DBethe^2-\[Epsilon]^2], 
	CompilationTarget->"C", RuntimeAttributes->{Listable}
];

(* dispersion relation for a d-dimensional hypercubic lattice *)
DispersionHypercubic = Compile[
	{{k,_Real,1}, {t,_Real}},
	-2.*t*Sum[Cos[ka], {ka, k}],
	CompilationTarget -> "C", RuntimeAttributes -> {Listable}
];

(* Local Green Function *)
LocalGreenFunction[DBethe_, \[CapitalSigma]_, EdMode_, zlist_, OptionsPattern[]] := Module[
	{Gloc, Floc, zero, d\[Epsilon], LocalGF, BrillouinZone, dk,
	Lattice = OptionValue[Lattice], d = OptionValue[LatticeDimension], LE = OptionValue[NumberOfPoints]},
	d\[Epsilon] = 2.*DBethe/LE;
	Which[
		EdMode == "Normal" && Lattice == "Bethe",
		Gloc = d\[Epsilon]*Total@Table[
			DoSBethe[\[Epsilon], DBethe]/(zlist - \[Epsilon] - \[CapitalSigma])
		, {\[Epsilon], -DBethe, DBethe, d\[Epsilon]}],
	(* -------------------------------------------------------- *)	
		EdMode == "Normal" && Lattice == "Hypercubic",
		LE = Floor[LE^(1/d)]; (* number of lattice points per size of the BZ *)
		dk = 2.Pi/LE; 
		BrillouinZone = Flatten[
			Outer[{##}&,##]&@@
				ConstantArray[
					Table[k, {k, -1.*Pi, 1.*Pi-dk, dk}],
					d
				], d-1];
		Gloc = (1./LE^d)*Total@Table[
			1./(zlist - DispersionHypercubic[k, DBethe] - \[CapitalSigma])
		, {k, BrillouinZone}],
	(* --------------------------------------------------- *)	
		EdMode == "Superc" && Lattice == "Bethe",
		Gloc = d\[Epsilon]*Total@Table[
				DoSBethe[\[Epsilon], DBethe]*(-zlist - Conjugate@\[CapitalSigma][[All,1,1]] - \[Epsilon])/(Abs[zlist - \[CapitalSigma][[All,1,1]] - \[Epsilon]]^2 + Abs[\[CapitalSigma][[All,1,2]]]^2)
			,{\[Epsilon], -DBethe, DBethe, d\[Epsilon]}];
		Floc = -\[CapitalSigma][[All,1,2]]*d\[Epsilon]*Total@Table[
				DoSBethe[\[Epsilon], DBethe]*(1./(Abs[zlist - \[CapitalSigma][[All,1,1]] - \[Epsilon]]^2 + Abs[\[CapitalSigma][[All,1,2]]]^2))
			,{\[Epsilon], -DBethe, DBethe, d\[Epsilon]}];
		LocalGF = Partition[#,2]&/@({Gloc, Floc, Conjugate@Floc, -Conjugate@Gloc}\[Transpose]),
	(* --------------------------------------------------- *)
		EdMode == "Superc" && Lattice == "Hypercubic",
		LE = Floor[LE^(1/d)]; (* number of lattice points per size of the BZ *)
		dk = 2.Pi/LE; 
		BrillouinZone = Flatten[
			Outer[{##}&,##]&@@
			ConstantArray[
				Table[k, {k, -1.*Pi, 1.*Pi-dk, dk}],
				d
			], d-1];
		Gloc = (1./LE^d)*Total@Table[
			-(zlist + Conjugate@\[CapitalSigma][[All,1,1]] + DispersionHypercubic[k, DBethe])/(Abs[zlist - \[CapitalSigma][[All,1,1]] - DispersionHypercubic[k, DBethe]]^2 + Abs[\[CapitalSigma][[All,1,2]]]^2)
		,{k, BrillouinZone}];
		Floc = -\[CapitalSigma][[All,1,2]]*(1./LE^d)*Total@Table[
			DispersionHypercubic[k, DBethe]*(1./(Abs[zlist - \[CapitalSigma][[All,1,1]] - DispersionHypercubic[k, DBethe]]^2 + Abs[\[CapitalSigma][[All,1,2]]]^2))
		,{k, BrillouinZone}];
		LocalGF = Partition[#,2]&/@({Gloc, Floc, Conjugate@Floc, -Conjugate@Gloc}\[Transpose])
	]
];
Options[LocalGreenFunction] = {Lattice -> "Bethe", LatticeDimension -> 2, NumberOfPoints -> 1000};



(* extract independent bath parameters depending on EdMode *)
TakeIndependentParameters[L_, f_, Norb_, \[Sigma]_, orb_, BathParameters_, EdMode_] :=
	Which[
		EdMode == "Normal",
		Join[
			BathParameters[[1]][[f*(orb-1)+\[Sigma]]], (* e1, e2, e3, ... *)
			BathParameters[[2]][[f*(orb-1)+\[Sigma]]] (* V1, V2, V3, ... *)
		],
	(* ---------------------------------------- *)
		EdMode == "Superc",
		Join[
			BathParameters[[1]][[f*(orb-1)+\[Sigma]]], (* e1, e2, e3, ... *)
			BathParameters[[2]][[f*(orb-1)+\[Sigma]]], (* V1, V2, V3, ... *)
			BathParameters[[3]][[f*(orb-1)+\[Sigma]]] (* \[CapitalDelta]1, \[CapitalDelta]2, \[CapitalDelta]3, ... *)
		],
	(* --------------------------------------- *)
		EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc",
		Flatten[BathParameters]
	];

(* correctly reshape the flat list of independent bath parameters *)
ReshapeBathParameters[L_, f_, Norb_, IndependentParameters_, EdMode_] := 
	Which[
		EdMode == "Normal",
		{ConstantArray[Take[IndependentParameters, L-1], f*Norb], (* e *)
		ConstantArray[Take[IndependentParameters, {L, 2(L-1)}], f*Norb]}, (* V *)
	(* ----------------------------------------- *)
		EdMode == "Superc",
		{ConstantArray[Take[IndependentParameters, L-1], f*Norb], (* e *)
		ConstantArray[Take[IndependentParameters, {L, 2(L-1)}], f*Norb], (* V *)
		ConstantArray[Take[IndependentParameters, {2L-1, 3(L-1)}], Norb]}, (* \[CapitalDelta] *)
	(* ----------------------------------------- *)
		EdMode == "InterorbNormal",
		{Partition[Take[IndependentParameters, (L-1)*f*Norb], L-1], (* e *)
		Partition[Take[IndependentParameters, {1+(L-1)*f*Norb, 2*(L-1)*f*Norb}], L-1]}, (* V *)
	(* ----------------------------------------- *)
		EdMode == "InterorbSuperc",
		{Partition[Take[IndependentParameters, (L-1)*f*Norb], L-1], (* e *)
		Partition[Take[IndependentParameters, {1+(L-1)*f*Norb, 2*(L-1)*f*Norb}], L-1], (* V *)
		Take[IndependentParameters, {1+2*(L-1)*f*Norb, 2*(L-1)*f*Norb + L-1}]}, (* \[CapitalXi] *)
	(* ----------------------------------------- *)
		EdMode == "FullSuperc",
		{Partition[Take[IndependentParameters, (L-1)*f*Norb], L-1], (* e *)
		Partition[Take[IndependentParameters, {1+(L-1)*f*Norb, 2*(L-1)*f*Norb}], L-1], (* V *)
		Partition[Take[IndependentParameters, {1+2*(L-1)*f*Norb, 1+2*(L-1)*f*Norb + (L-1)*Norb}], L-1], (* \[CapitalDelta] *)
		IndependentParameters[[1+2*(L-1)*f*Norb + (L-1)*Norb;;]]} (* \[CapitalXi] *)
	];

(* Perform minimization according to the self consistency condition *)
SelfConsistency[DBethe_, \[Mu]_, Weiss_, symbols_, z_, IndependentParameters_, LocalG_, zlist_, EdMode_, OptionsPattern[]] := Module[
	{Lattice = OptionValue[Lattice],
	LFit = Min[Length[zlist], OptionValue[NumberOfFrequencies]],
	residue, newparameters, \[Chi]},
	(* define the target function to minimize depending on the Lattice and EdMode (if Lattice = Bethe there is a shortcut) *)
	Which[
		Lattice == "Bethe",
		\[Chi][symbols] = Mean[Abs[
			((Weiss - z - \[Mu])/.{z -> Take[zlist, LFit]}) + (DBethe^2/4.)*Take[LocalG, LFit]
		]^2]
	];
	{residue, newparameters} =
		FindMinimum[
			\[Chi][symbols],
			{symbols, StartingParameters}\[Transpose],
			Method -> "ConjugateGradient",
			MaxIterations -> 700,
			AccuracyGoal -> 5
		];
	symbols/.newparameters
];
Options[SelfConsistency] = {Lattice -> "Bethe", NumberOfFrequencies -> 2000};



End[];

EndPackage[];
