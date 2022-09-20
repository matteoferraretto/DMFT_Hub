(* ::Package:: *)

(* ::Subtitle:: *)
(*Official functions*)


(* check if pair hopping is possible *)
PairHoppingQ = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer},{orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	If[
		(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+1]])[[i]]==0&&
	(IntegerDigits[#,2,L]&@state[[f*(orb1-1)+2]])[[i]]==0&&
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+1]])[[i]]==1&&
	(IntegerDigits[#,2,L]&@state[[f*(orb2-1)+2]])[[i]]==1,
		True,
	(*else*)
		False
	]
];

(* selects states for which pair creation (i,orb1,\[Sigma]1), (j,orb2,\[Sigma]2) is possible*)
PairHoppingSelect = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {stateList,_Integer,2}
	},
	Select[stateList, PairHoppingQ[L,f,i,orb1,orb2,#]&],
	RuntimeAttributes->{Listable}, Parallelization->True, CompilationTarget->"C"
];

(* pair hopping of a state *)
PairHopping = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	MapAt[
		BitOr[#,2^(L-i)]&,
		MapAt[
			BitOr[#,2^(L-i)]&,
			MapAt[
				BitAnd[#, BitNot[-2^(L-i)]]&,
				MapAt[
					BitAnd[#, BitNot[-2^(L-i)]]&,
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



(* ::Subtitle:: *)
(*Example / Test*)


L=2; f=2; Norb=2;

(* take two states *)
\[Psi]={3,3,0,0};
\[Xi]={0,0,2,2};
IntegerDigits[#,2,L]&@\[Psi]
IntegerDigits[#,2,L]&@\[Xi]

(* check if pair hopping is possible *)
PairHoppingQ[L,f,1,2,1,\[Psi]]
PairHoppingQ[L,f,1,2,1,\[Xi]]

(* select those states for which pair hopping is possible *)
PairHoppingSelect[L,f,1,2,1,#]&@{\[Psi],\[Xi]}

(* when possible, apply pair hopping operator *)
IntegerDigits[#,2,L]&@PairHopping[L,f,1,2,1,\[Psi]]


(* ::Subtitle:: *)
(**)
(*Performance*)


L=2; f=2; Norb=2;

\[Psi]={3,3,0,0};

AbsoluteTiming[
	PairHopping[L,f,1,2,1,#]&/@ConstantArray[\[Psi],300000];
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

(* pair hopping *)
PairHoppingAlternative = Compile[{
	{L,_Integer}, {f,_Integer}, {i,_Integer}, {orb1,_Integer}, {orb2,_Integer}, {state,_Integer,1}
	},
	cdg[L,f,i,1,orb1,#]&@cdg[L,f,i,2,orb1,#]&@c[L,f,i,2,orb2,#]&@c[L,f,i,1,orb2,#]&@state,
	CompilationTarget->"C"
];

(* check performance and see it's worse *)
L=2; f=2; Norb=2;
\[Psi]={3,3,0,0};

AbsoluteTiming[
	PairHoppingAlternative[L,f,1,2,1,#]&/@ConstantArray[\[Psi],300000];
]
