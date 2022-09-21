(* ::Package:: *)

(* ::Subtitle:: *)
(*Official functions*)


(* Counts how many fermions there are before state[[\[Sigma],orb,i]] (excluded). If you set i=1,\[Sigma]=1,orb=1 you get 0; if you set i=L+1,\[Sigma]=f,orb=Norb you get the total number of fermions in the state *)
CountFermions = Compile[{
	{L,_Integer},{f,_Integer},{i,_Integer},{\[Sigma],_Integer},{orb,_Integer},{state,_Integer,1}
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



(* ::Subtitle:: *)
(*Examples / Test*)


L=2; f=2; Norb=2;
i=1; \[Sigma]=1; orb=2;
j=2;

(* define two states *)
\[Psi]={1,2,2,1};
\[Xi]={3,1,1,1};
IntegerDigits[#,2,L]&@\[Psi]
IntegerDigits[#,2,L]&@\[Xi]

(* count how many fermions there are before position (i,\[Sigma],orb) *)
CountFermions[L,f,i,\[Sigma],orb,\[Psi]]
CountFermions[L,f,i,\[Sigma],orb,\[Xi]]

(* show the sign obtained moving a fermionic operator in front of the position i,\[Sigma],orb *)
CSign[L,f,i,\[Sigma],orb,\[Psi]]
CSign[L,f,i,\[Sigma],orb,\[Xi]]

(* show the sign obtained moving a fermionic operator in front of the position i,\[Sigma],orb and another in front of position j,\[Sigma],orb *)
CCSign[L,f,{i,j},{\[Sigma],\[Sigma]},{orb,orb},\[Psi]]
CCSign[L,f,{i,j},{\[Sigma],\[Sigma]},{orb,orb},\[Xi]]


(* ::Subtitle:: *)
(*Performance*)


L=2; f=2; Norb=2;
i=1; \[Sigma]=1; orb=2;
\[Psi]={1,2,2,1};

AbsoluteTiming[
	CountFermions[L,f,i,\[Sigma],orb,#]&/@ConstantArray[\[Psi],300000];
]
AbsoluteTiming[
	CSign[L,f,i,\[Sigma],orb,#]&/@ConstantArray[\[Psi],300000];
]
AbsoluteTiming[
	CCSign[L,f,{i,j},{\[Sigma],\[Sigma]},{orb,orb},#]&/@ConstantArray[\[Psi],300000];
]


(* ::Subtitle:: *)
(*Alternative formulations*)


CountFermionsAlternative[L_,f_,i_,\[Sigma]_,orb_,state_]:=
	Total@Take[ (* sum *)
		Flatten[IntegerDigits[#,2,L]&/@state] (* flattened version of the binary representation of the state *)
		,L*(f*(orb-1)+\[Sigma]-1)+(i-1) (* sum up to this index in the flattened version of the state *)
	];
CountFermionsAlternative[L_,f_,i_,j_,\[Sigma]1_,\[Sigma]2_,orb1_,orb2_,state_]:=
	Total@Take[ (* sum *)
		Flatten[IntegerDigits[#,2,L]&/@state] (* flattened version of the binary representation of the state *)
		,{L*(f*(orb1-1)+\[Sigma]1-1)+i, L*(f*(orb2-1)+\[Sigma]2-1)+(j-1)} (* sum up to this index in the flattened version of the state *)
	];
	
CSignAlternative[L_, f_, i_, \[Sigma]_, orb_, state_]:=(-1)^CountFermions[L,f,i,\[Sigma],orb,state];

CCSignAlternative[L_, f_, i_, j_, \[Sigma]1_, \[Sigma]2_, orb1_, orb2_, state_] := (-1)^CountFermions[L, f, i, j, \[Sigma]1, \[Sigma]2, orb1, orb2, state]
(* CCSign just works with 2 operators ... *)

L=2; f=2; Norb=2;
i=1; \[Sigma]=1; orb=2;
\[Psi]={1,2,2,1};
	
AbsoluteTiming[
	CountFermionsAlternative[L,f,i,\[Sigma],orb,#]&/@ConstantArray[\[Psi],300000];
]

AbsoluteTiming[
	CSignAlternative[L,f,i,\[Sigma],orb,#]&/@ConstantArray[\[Psi],300000];
]



