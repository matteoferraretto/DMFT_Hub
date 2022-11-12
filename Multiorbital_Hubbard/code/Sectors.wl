(* ::Package:: *)

BeginPackage["Sectors`"]


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


Begin["Private`"]

Print["Package Sectors` loaded successfully."];

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

End[]

EndPackage[]
