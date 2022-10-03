(* ::Package:: *)

BeginPackage["DMFT`"]


StartingBath::usage = "StartingBath[InitializeBathMode, Nbath, Norb, EdMode] returns a list containing the bath parameters to start the DMFT loop.
If EdMode = ''Normal'' then the output has the form {e,V}, where e and V are lists of Norb x Nbath elements, representing the bath energies and the bath-impurity hybridizations.
If EdMode = ''Superc'' then the output has the form {e,V,\[CapitalDelta]}, where e and V are defined as above, and \[CapitalDelta] is the Norb x Nbath dimensional list of pairs creation (annihilation) amplitudes.
If EdMode = ''InterorbSuperc'' then the output has the form {e,V,\[CapitalDelta],\[CapitalXi]}, where e, V, \[CapitalDelta] are as above, and \[CapitalXi] is the Nbath - dimensional list of interorbital pairs creation (annihilation) amplitudes
InitializeBathMode is a string with the path to the file containing the bath parameters; if it is set to ''Default'', default parameters are dropped."


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


n::usage = "n[L, f, Norb, i, \[Sigma], orb, state]"


(* Hamiltonian defining functions *)
HNonlocal::usage = "HNonlocal[L, f, Norb, Sectors, EdMode]"

HNonlocalInfo::usage = "HNonlocalInfo[L, f, Norb, EdMode] prints useful info about the order of non local Hamiltonian blocks."

HLocal::usage = "HLocal[L, f, Norb, Sectors, EdMode] "

HImp::usage = "..."


Begin["`Private`"];


(* Initialize starting bath *)
StartingBath[InitializeBathMode_, Nbath_, Norb_, EdMode_] := Module[
	{e,V,\[CapitalDelta],\[CapitalXi]},
	Which[
		EdMode == "Normal" || EdMode == "InterorbNormal",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb],	
		(*else*)
			{e,V} = Import[InitializeBathMode,"Table"];
		];
		Return[{e,V}],
(* ---------------------------------------------- *)
		EdMode == "Superc",
		If[
			InitializeBathMode == "Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb];
			\[CapitalDelta] = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb],
		(*else*)
			{e,V,\[CapitalDelta]} = Import[InitializeBathMode,"Table"];
		];
		Return[{e,V,\[CapitalDelta]}],
(* ---------------------------------------------- *)
		EdMode == "InterorbSuperc",
		If[
			InitializeBathMode=="Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb];
	        \[CapitalXi] = Table[1.,{k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalXi]} = Import[InitializeBathMode,"Table"];
		];
	   Return[{e,V,\[CapitalXi]}],
(* ---------------------------------------------- *)
		EdMode == "FullSuperc",
		If[
			InitializeBathMode=="Default",
			e = ConstantArray[
				Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
			Norb];
			V = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb];
			\[CapitalDelta] = ConstantArray[
				Table[1.,{k,1,Nbath}],
			Norb];
	        \[CapitalXi] = Table[1.,{k,1,Nbath}],
		(*else*)
			{e,V,\[CapitalDelta],\[CapitalXi]} = Import[InitializeBathMode,"Table"];
		];
	   Return[{e,V,\[CapitalDelta],\[CapitalXi]}];
	]
];

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
BASIS[L_Integer, flavors_Integer, qns_]:=Flatten[Outer[{##}&,##]&@@Table[basis[L,qns[[\[Sigma]]]],{\[Sigma],1,flavors}],flavors-1];

(* list of all the quantum numbers that label the sectors *)
SectorList[L_,f_,Norb_,EdMode_]:=Module[
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
DimSector[L_, f_, Norb_, qns_, EdMode_String]:=Module[
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
		,{ndw,Max[-qns,0],Min[Norb*L-qns,Norb*L]}];
	];
	dim
];

(* build sector, i.e. list of all the Fock states with given quantum number(s) qns *)
BuildSector[L_, f_, Norb_, qns_, EdMode_]:=Module[
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
DrawState[L_, f_, Norb_, j_, \[Sigma]_, orb_]:=Module[
	{impuritycolor, color},
	impuritycolor[i_]:=If[i==1,Blue,Black];
	color[i_,index_]:=If[i==j && index==f*(orb-1)+\[Sigma], Red, Black];
	Table[
		Graphics[
			Table[
				{EdgeForm[{Thick,impuritycolor[i]}],color[i,index],Opacity[.3],Rectangle[{1.1*i,0}]},{i,L}
			]
		],
	{index,f*Norb}]
];
DrawState[L_,f_,Norb_]:=Module[{},
	Print[Style["Architecture of a state",16]];
	Print[Style["Red:",Red]," (site j, spin \[Sigma], orbital orb)"];
	Print[Style["Blue edge:",Blue]," impurity"];
	Manipulate[
		DrawState[L,f,Norb,j,\[Sigma],orb],
	{j,1,L,1}, {\[Sigma],1,f,1}, {orb,1,Norb,1}]
];


(*               MANIPULATION OF STATES                *)
(* apply cdg_orb_spin in the impurity site to a basis state: return the integer form of the resulting basis state or 0 *)
cdg = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}
	},
	MapAt[
		BitSet[#, L-i]&,
		state,
		f*(orb-1)+\[Sigma]
	], CompilationTarget->"C"];

(* apply c_orb_spin in the impurity site to a basis state: return the integer form of the resulting basis state or 0 *)
c = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}
	},
	MapAt[
		BitClear[#, L-i]&,
		state,
		f*(orb-1)+\[Sigma]
	], CompilationTarget->"C"];

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
	(-1)^CountFermions[L,f,i,\[Sigma],orb,state],
	CompilationTarget->"C"
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
	], CompilationTarget->"C"
];

CCCCSign = CCSign;


(*               HOPPING FUNCTIONS             *)
(* gives True if hopping from (j, \[Sigma]2, orb2) to (i,\[Sigma]1,orb1) is possible, False otherwise *)
HopQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {j,_Integer}, {\[Sigma]1,_Integer}, {\[Sigma]2,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	If[(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+\[Sigma]2]])[[j]]==1 && (IntegerDigits[#,2,L]&@state[[f*(orb1-1)+\[Sigma]1]])[[i]]==0,
		True,
	(* else *)
		False
	], 
	RuntimeAttributes->{Listable}, Parallelization->True
];

(* select states for which hopping (j, \[Sigma]2, orb2) \[Rule] (i, \[Sigma]1, orb1) is possible *)
HopSelect = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {j,_Integer}, {\[Sigma]1,_Integer}, {\[Sigma]2,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {stateList,_Integer,2}
	},
	Select[stateList, HopQ[L,f,i,j,\[Sigma]1,\[Sigma]2,orb1,orb2,#]&],
	RuntimeAttributes->{Listable}, Parallelization->True, CompilationTarget->"C"
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
	CompilationTarget->"C"
];

(* number of particles on site i with spin \[Sigma] in orbital orb *)
n[L_, f_, Norb_, i_, \[Sigma]_, orb_, state_] := (IntegerDigits[#,2,L]&@state[[f*(orb-1)+\[Sigma]]])[[i]];


(*           PAIR CREATION / ANNIHILATION FUNCTIONS          *)
(* gives True if it is possible to create a pair (i,orb1,\[Sigma]1) and (j,orb2,\[Sigma]2) and False otherwise *)
CreatePairQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {j,_Integer}, {\[Sigma]1,_Integer}, {\[Sigma]2,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	If[
		(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+\[Sigma]1]])[[i]]==0&&(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+\[Sigma]2]])[[j]]==0,
		True,
	(*else*)
		False
	]
];

(* selects states for which pair creation (i,orb1,\[Sigma]1), (j,orb2,\[Sigma]2) is possible*)
CreatePairSelect = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {j,_Integer}, {\[Sigma]1,_Integer}, {\[Sigma]2,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {stateList,_Integer,2}
	},
	Select[stateList, CreatePairQ[L,f,i,j,\[Sigma]1,\[Sigma]2,orb1,orb2,#]&],
	RuntimeAttributes->{Listable}, Parallelization->True, CompilationTarget->"C"
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
	CompilationTarget->"C"
];


(*          PAIR HOPPING FUNCTIONS         *)
(* gives True if it is possible to make a pair hopping from (i,orb2) to (i,orb1), False otherwise  *)
PairHoppingQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer},{orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	If[
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+1]])[[i]] == 0 &&
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+2]])[[i]] == 0 &&
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+1]])[[i]] == 1 &&
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+2]])[[i]] == 1,
		True,
	(*else*)
		False
	]
];

(* selects states for which pair hopping from (i,orb2) to (i,orb1) is possible*)
PairHoppingSelect = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {stateList,_Integer,2}
	},
	Select[stateList, PairHoppingQ[L,f,i,orb1,orb2,#]&],
	RuntimeAttributes->{Listable}, Parallelization->True, CompilationTarget->"C"
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
	CompilationTarget->"C"
];


(*           SPIN EXCHANGE FUNCTIONS         *)
(* gives True if spin exchange on site i is possible, False otherwise: spin exchange = cdg_orb1_up c_orb1_dw cdg_orb2_dw c_orb2_up *)
SpinExchangeQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	If[
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+1]])[[i]] == 1 && 
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+2]])[[i]] == 1 &&
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+2]])[[i]] == 0 && 
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+1]])[[i]] == 0,
		True,
	(* else *)
		False
	], 
	RuntimeAttributes->{Listable}, Parallelization->True
];

(* select states for which spin exchange on site i is possible: spin exchange = cdg_orb1_up c_orb1_dw cdg_orb2_dw c_orb2_up *)
SpinExchangeSelect = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {stateList,_Integer,2}
	},
	Select[stateList, SpinExchangeQ[L,f,i,orb1,orb2,#]&],
	RuntimeAttributes->{Listable}, Parallelization->True, CompilationTarget->"C"
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
	CompilationTarget->"C"
];


(*          BUILD THE IMPURITY HAMILTONIAN           *)
(* Non-local Hamiltonian blocks *)
HNonlocal[L_, f_, Norb_, Sectors_, EdMode_] := Module[
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
				Do[
					num = n[L,f,Norb,j,\[Sigma],orb,#]&/@\[Psi];(*local density*)
					Hblock += SparseArray@DiagonalMatrix[num];
				,{\[Sigma],1,f}];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Hopping",
				Do[
					\[Psi]1 = HopSelect[L, f, 1, j, \[Sigma], \[Sigma], orb, orb, \[Psi]];
					If[Length[\[Psi]1] == 0, Continue[];];
					\[Chi] = Hop[L,f,1,j,\[Sigma],\[Sigma],orb,orb,#]&/@(\[Psi]1);
					rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
					\[CapitalSigma] = (CCSign[L,f,{1,j},{\[Sigma],\[Sigma]},{orb,orb},#]&/@\[Psi]1);
					Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
				,{\[Sigma],1,f}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Superc" && (EdMode == "Superc" || EdMode == "FullSuperc"),
				\[Psi]1 = CreatePairSelect[L,f,j,j,1,2,orb,orb,#]&@\[Psi];
				\[Chi] = CreatePair[L,f,j,j,1,2,orb,orb,#]&/@\[Psi]1;
				rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
				\[CapitalSigma] = (CCSign[L,f,{j,j},{1,2},{orb,orb},#]&/@\[Psi]1);
				Hblock += SparseArray[pos->\[CapitalSigma], {dim,dim}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "InterorbSuperc" && (EdMode == "InterorbSuperc" || EdMode == "FullSuperc"),
				Do[
					If[
						orb2 > orb,
						\[Psi]1 = CreatePairSelect[L,f,j,j,1,2,orb,orb2,#]&@\[Psi];
						\[Chi] = CreatePair[L,f,j,j,1,2,orb,orb2,#]&/@\[Psi]1;
						rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						\[CapitalSigma] = (CCSign[L,f,{j,j},{1,2},{orb,orb2},#]&/@\[Psi]1);
						Hblock += SparseArray[pos->\[CapitalSigma], {dim,dim}];
					]
				, {orb2,1,Norb}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];
			];
		,{flag, {"Bath","Hopping","Superc","InterorbSuperc"}}, {orb,1,Norb}, {j,2,L}];
		AppendTo[H, Hsector];
	,{\[Psi], Sectors}];
	H
];

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
HLocal[L_, f_, Norb_, Sectors_, EdMode_] := Module[
	{H,Hsector,Hblock,rules,dispatch,num,dim,\[Psi]1,\[Chi],rows,cols,pos,\[CapitalSigma],Nimp=1},
	H = {};
	Do[
		Hsector = {};
		dim = Length[\[Psi]];
		Hsector = {};
		rules = Flatten[MapIndexed[{#1->#2[[1]]}&,\[Psi]],1];
		dispatch = Dispatch[rules];
		Do[
			Hblock = SparseArray[{},{dim,dim}];
			Which[
				flag == "Hubbard",
				Do[
					num = Sum[
						n[L,f,Norb,j,1,orb,#]*n[L,f,Norb,j,2,orb,#]
					,{j,1,Nimp}]&/@\[Psi];
					Hblock = SparseArray@DiagonalMatrix[num];
					AppendTo[Hsector, Hblock];
				,{orb,1,Norb}],
			(* ---------------------------------- *)
				flag == "Interorb_Hubbard_Opposite_Spin" && Norb > 1,
				num = Sum[
					If[orbA != orbB,
						n[L,f,Norb,j,1,orbA,#]*n[L,f,Norb,j,2,orbB,#],
					(*else*)
						0]
				,{orbA,1,Norb}, {orbB,1,Norb}, {j,1,Nimp}]&/@\[Psi];
				Hblock = SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];,
			(* ---------------------------------- *)
				flag == "Interorb_Hubbard_Same_Spin" && Norb > 1,
				num = Sum[
					n[L,f,Norb,j,\[Sigma],1,#]*n[L,f,Norb,j,\[Sigma],2,#]
				,{\[Sigma],1,f}, {j,1,Nimp}]&/@\[Psi];
				Hblock = SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector,Hblock];,
			(* ---------------------------------- *)
				flag == "Pair_Hopping" && Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "Superc" || EdMode == "FullSuperc"),
				Do[
					If[orb2 > orb1,
						\[Psi]1 = PairHoppingSelect[L, f, 1, orb1, orb2, #]&@\[Psi];
						If[Length[\[Psi]1] == 0, Continue[];];
						\[Chi] = PairHopping[L, f, 1, orb1, orb2, #]&/@\[Psi]1;
						rows=\[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
						\[CapitalSigma]=CCSign[L, f, {j,j,j,j}, {1,2,1,2}, {orb1,orb1,orb2,orb2}, #]&/@\[Psi]1;
						Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					];
				,{orb1,1,Norb}, {orb2,1,Norb}, {j,1,Nimp}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector,Hblock];,
			(* ----------------------------------- *)
				flag == "Spin_Exchange" && Norb > 1 && (EdMode == "InterorbNormal" || EdMode == "InterorbSuperc" || EdMode == "FullSuperc"),
				Do[
					If[orb2 > orb1,
						\[Psi]1 = SpinExchangeSelect[L, f, 1, orb1, orb2, #]&@\[Psi];
						If[Length[\[Psi]1] == 0, Continue[];];
						\[Chi] = SpinExchange[L, f, 1, orb1, orb2, #]&/@\[Psi]1;
						rows = \[Chi]/.dispatch;(* *)cols = \[Psi]1/.dispatch;(* *)pos = {rows,cols}\[Transpose];
						\[CapitalSigma] = CCSign[L, f, {j,j,j,j}, {1,2,1,2}, {orb1,orb1,orb2,orb2}, #]&/@\[Psi]1;
						Hblock += SparseArray[pos->\[CapitalSigma],{dim,dim}];
					];
				,{orb1,1,Norb}, {orb2,1,Norb}, {j,1,Nimp}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];,
			(* ---------------------------------- *)
				flag == "Energy_Shift",
				num = Sum[
					n[L,f,Norb,j,\[Sigma],orb,#]
				,{\[Sigma],1,f}, {orb,1,Norb}, {j,1,Nimp}]&/@\[Psi];
				Hblock = SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];
			];
		,{flag, {"Hubbard","Interorb_Hubbard_Opposite_Spin","Interorb_Hubbard_Same_Spin","Pair_Hopping","Spin_Exchange","Energy_Shift"}}];
		AppendTo[H,Hsector];
	,{\[Psi],Sectors}];
	H
];


HImp[L_, f_, Norb_, Sectors_, BathParameters_, InteractionParameters_, EdMode_] := Module[{
	HnonlocBlocks = HNonlocal[L, f, Norb, Sectors, EdMode],
	HlocBlocks = HLocal[L, f, Norb, Sectors, EdMode],
	EffectiveInteractionParameters,
	Hloc, Hnonloc
	},
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
			BathParameters[[i]]*#[[i]],
		{i, 1, Length@BathParameters}]&/@HnonlocBlocks
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
(* gives right neighbor of a site j in a closed (open) chain of L elements *)
Neighbor = Compile[{
	{L,_Integer}, {j,_Integer}
	},
	Mod[j,L]+1
];

Options[Hnonint] = {RealPBC -> True, SyntheticPBC -> False, RealPhase -> {0}, SyntheticPhase -> 0};

Hnonint[L_, f_, Norb_:1, Sectors_, EdMode_, OptionsPattern[]] := Module[
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
				num = n[L, f, Norb, j, \[Sigma], orb, #]&/@\[Psi];(*local density*)
				Hblock += SparseArray@DiagonalMatrix[num];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Hopping",
				If[j == L && !OptionValue[RealPBC], Continue[];];
				\[Psi]1 = HopSelect[L, f, j, Neighbor[L,j], \[Sigma], \[Sigma], orb, orb, \[Psi]];
				If[Length[\[Psi]1] != 0, Continue[];];
				\[Chi] = Hop[L, f, j, Neighbor[L,j], \[Sigma], \[Sigma], orb, orb, #]&/@(\[Psi]1);
				rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
				\[CapitalSigma] = (CCSign[L,f,{j,Neighbor[L,j]},{\[Sigma],\[Sigma]},{orb,orb},#]&/@\[Psi]1);
				If[
					OptionValue[RealPhase] == {0},
					Hblock += SparseArray[pos -> \[CapitalSigma],{dim,dim}];,
				(* else *)
					Hblock += SparseArray[pos -> \[CapitalSigma]*Exp[I*OptionValue[RealPhase][[\[Sigma]]]],{dim,dim}]
					]
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];,
			(* --------------------------------------------------------------- *)
				flag == "Raman",
				If[\[Sigma] == f && !OptionValue[SyntheticPBC], Continue[];];
				\[Psi]1 = HopSelect[L, f, j, j, \[Sigma], Neighbor[f,\[Sigma]], orb, orb, \[Psi]];
				If[Length[\[Psi]1] != 0, Continue[];];
				\[Chi] = Hop[L, f, j, j, \[Sigma], Neighbor[f,\[Sigma]], orb, orb, #]&/@(\[Psi]1);
				rows = \[Chi]/.dispatch;(* *)cols=\[Psi]1/.dispatch;(* *)pos={rows,cols}\[Transpose];
				\[CapitalSigma] = (CCSign[L,f,{j,j},{\[Sigma],Neighbor[f,\[Sigma]]},{orb,orb},#]&/@\[Psi]1);
				Hblock += SparseArray[pos -> \[CapitalSigma]*Exp[I*OptionValue[SyntheticPhase]], {dim,dim}];
				Hblock = Hblock + Hblock\[ConjugateTranspose];
				AppendTo[Hsector, Hblock];	
			];
		,{flag, {"Potential","Hopping","Raman"}}, {orb,1,Norb}, {\[Sigma],1,f}, {j,1,L}];
		AppendTo[H, Hsector];
	,{\[Psi], Sectors}];
	H
];


End[];

EndPackage[];
