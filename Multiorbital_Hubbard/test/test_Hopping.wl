(* ::Package:: *)

(* ::Subtitle:: *)
(*Official functions*)


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

(* select states for which hopping (j, \[Sigma]1, orb1) \[Rule] (i, \[Sigma]2, orb2) is possible *)
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
		BitOr[#,2^(L-i)]&,
		MapAt[
			BitAnd[#, BitNot[-2^(L-j)]]&,
			state,
			f*(orb1-1)+\[Sigma]1
		],
		f*(orb2-1)+\[Sigma]2
	],
	CompilationTarget->"C"
];



(* ::Subtitle:: *)
(*Examples / Test*)


L=2; f=2; Norb=2;
i=2; j=1; \[Sigma]=2; orb=2;

(* define two states *)
\[Psi]={1,1,2,2};
\[Xi]={0,0,1,1};
IntegerDigits[#,2,L]&@\[Psi]
IntegerDigits[#,2,L]&@\[Xi]

(* check if hopping is possible *)
HopQ[L,f,i,j,\[Sigma],\[Sigma],orb,orb,#]&@\[Psi]
HopQ[L,f,i,j,\[Sigma],\[Sigma],orb,orb,#]&@\[Xi]

(* select the states where the given hopping process is possible *)
{\[Psi],\[Xi]}
HopSelect[L,f,i,j,\[Sigma],\[Sigma],orb,orb,#]&@{\[Psi],\[Xi]}

(* apply hopping operator to the state *)
IntegerDigits[#,2,L]&@\[Psi]
IntegerDigits[#,2,L]&@(Hop[L,f,i,j,\[Sigma],\[Sigma],orb,orb,\[Psi]])


(* ::Subtitle:: *)
(**)
(*Performance*)


L=2; f=2; Norb=2;
i=2; j=1; \[Sigma]=2; orb=2;
\[Psi]={1,1,2,2};

AbsoluteTiming[
	Hop[L,f,i,j,\[Sigma],\[Sigma],orb,orb,#]&/@ConstantArray[\[Psi],300000];
]


(* ::Subtitle:: *)
(**)
(*Alternative formulations*)


cdg = Compile[{
	{L,_Integer}, {f,_Integer}, {j,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}
	},
	MapAt[
		BitOr[#, 2^(L-j)]&,
		state,
		f*(orb-1)+\[Sigma]
	], CompilationTarget->"C"];

(* apply c_orb_spin in the impurity site to a basis state: return the integer form of the resulting basis state or 0 *)
c = Compile[{
	{L,_Integer}, {f,_Integer}, {j,_Integer}, {\[Sigma],_Integer}, {orb,_Integer}, {state,_Integer,1}
	},
	MapAt[
		BitAnd[#, BitNot[-2^(L-j)]]&,
		state,
		f*(orb-1)+\[Sigma]
	], CompilationTarget->"C"];
	
HopAlternative = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {j,_Integer}, {\[Sigma]1,_Integer}, {\[Sigma]2,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}},
	cdg[L,f,i,\[Sigma]1,orb1,#]&@c[L,f,j,\[Sigma]2,orb2,#]&@state
	, CompilationTarget->"C"];

IntegerDigits[#,2,L]&@\[Psi]
IntegerDigits[#,2,L]&@HopAlternative[L,f,i,j,\[Sigma],\[Sigma],orb,orb,\[Psi]]

AbsoluteTiming[
	HopAlternative[L,f,i,j,\[Sigma],\[Sigma],orb,orb,#]&/@ConstantArray[\[Psi],300000];
]



