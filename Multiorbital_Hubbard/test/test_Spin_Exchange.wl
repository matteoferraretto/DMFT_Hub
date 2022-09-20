(* ::Package:: *)

(* ::Subtitle:: *)
(*Official functions*)


(* gives True if hopping from (j, \[Sigma]2, orb2) to (i,\[Sigma]1,orb1) is possible, False otherwise *)
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

(* select states for which hopping (j, \[Sigma]1, orb1) \[Rule] (i, \[Sigma]2, orb2) is possible *)
SpinExchangeSelect = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {stateList,_Integer,2}
	},
	Select[stateList, SpinExchangeQ[L,f,i,orb1,orb2,#]&],
	RuntimeAttributes->{Listable}, Parallelization->True, CompilationTarget->"C"
];

(* apply hopping operator to a given state *)
SpinExchange = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}},
	MapAt[
		BitOr[#,2^(L-i)]&,
		MapAt[
			BitOr[#,2^(L-i)]&,
			MapAt[
				BitAnd[#, BitNot[-2^(L-i)]]&,
				MapAt[
					BitAnd[#, BitNot[-2^(L-i)]]&,
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



(* ::Subtitle:: *)
(*Examples / Test*)


L=2; f=2; Norb=2;
i=1;

(* define two states *)
\[Psi]={1,2,2,1};
\[Xi]={0,0,1,1};
IntegerDigits[#,2,L]&@\[Psi]
IntegerDigits[#,2,L]&@\[Xi]

(* check if hopping is possible *)
SpinExchangeQ[L,f,i,1,2,#]&@\[Psi]
SpinExchangeQ[L,f,i,1,2,#]&@\[Xi]

(* select the states where the given hopping process is possible *)
{\[Psi],\[Xi]}
SpinExchangeSelect[L,f,i,1,2,#]&@{\[Psi],\[Xi]}

(* apply hopping operator to the state *)
IntegerDigits[#,2,L]&@\[Psi]
IntegerDigits[#,2,L]&@(SpinExchange[L,f,i,1,2,\[Psi]])


(* ::Subtitle:: *)
(*Performance*)


L=2; f=2; Norb=2;
i=1;
\[Psi]={1,2,2,1};

AbsoluteTiming[
	SpinExchange[L,f,i,1,2,#]&/@ConstantArray[\[Psi],300000];
]


(* ::Subtitle:: *)
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

SpinExchangeAlternative = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}},
	cdg[L,f,i,1,orb1,#]&@c[L,f,i,2,orb1,#]&@cdg[L,f,i,2,orb2,#]&@c[L,f,i,1,orb2,#]&@state
	, CompilationTarget->"C"];

\[Psi]={1,2,2,1};
IntegerDigits[#,2,L]&@\[Psi]
IntegerDigits[#,2,L]&@SpinExchangeAlternative[L,f,i,1,2,\[Psi]]

AbsoluteTiming[
	SpinExchangeAlternative[L,f,i,1,2,#]&/@ConstantArray[\[Psi],300000];
]




