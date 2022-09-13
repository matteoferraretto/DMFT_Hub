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
If EdMode = ''Superc'' each element in the list has the form {sz_orb1, sz_orb2, ...} where sz_orb = n_orb_up-n_orb_dw (it only makes sense for spin 1/2 particles so far).
If EdMode = ''InterorbSuperc'' each element in the list has the form sz_tot where sz_tot = \!\(\*SubscriptBox[\(\[Sum]\), \(orb\)]\)(n_orb_up - n_orb_dw) (it only makes sense for spin 1/2 particles so far)."

BuildSector::usage = "BuildSector[L, f, qns, EdMode] creates a list of basis of the Fock subpace with given quantum numbers qns. 
If EdMode = ''Normal'', then qns={n,nup} are the particle number and the spin up particle number. If nup is set to Null, the function returns all the states with n particles and any value of nup.
If EdMode = ''Superc'', then qns=sz is the total spin-z of the fermions, i.e. nup - ndw.
The states are created in the integer represenation. L is the total number of sites, f is the number of flavours (only supports f=2 at the moment). "


Begin["`Private`"];


(* Initialize starting bath *)
StartingBath[InitializeBathMode_String, Nbath_Integer, Norb_Integer, EdMode_String]:=Module[
	{e,V,\[CapitalDelta],\[CapitalXi]},
	Which[
		EdMode=="Normal",
		If[
			InitializeBathMode=="Default",
			e=ConstantArray[
			Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
		Norb];
			V=ConstantArray[
			Table[1.,{k,1,Nbath}],
		Norb],	
		(*else*)
			{e,V}=Import[InitializeBathMode,"Table"];
		];
		Return[{e,V}],
	(* ---------------------------------------------- *)
		EdMode=="Superc",
		If[
			InitializeBathMode=="Default",
			e=ConstantArray[
			Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
		Norb];
			V=ConstantArray[
			Table[1.,{k,1,Nbath}],
		Norb];
			\[CapitalDelta]=ConstantArray[
			Table[1.,{k,1,Nbath}],
		Norb],
		(*else*)
	     {e,V,\[CapitalDelta]}=Import[InitializeBathMode,"Table"];
		];
	Return[{e,V,\[CapitalDelta]}],
(* ---------------------------------------------- *)
	EdMode=="InterorbSuperc",
	If[
		InitializeBathMode=="Default",
		e=ConstantArray[
			Table[-(Nbath-1)/2.+k,{k,0,Nbath-1}],
		Norb];
			   V=ConstantArray[
			Table[1.,{k,1,Nbath}],
		Norb];
			   \[CapitalDelta]=ConstantArray[
			Table[1.,{k,1,Nbath}],
		Norb];
	         \[CapitalXi]=Table[1.,{k,1,Nbath}],
	(*else*)
		{e,V,\[CapitalDelta],\[CapitalXi]}=Import[InitializeBathMode,"Table"];
	];
	   Return[{e,V,\[CapitalDelta],\[CapitalXi]}];
	]
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
SectorList[L_, f_, Norb_, EdMode_]:=Module[
	{QnsSectorList},
	Which[
		EdMode == "Normal",
		QnsSectorList=Flatten[
			Outer[{##}&,##]&@@ConstantArray[
			Range[0,L]
		,Norb*f],Norb*f-1],
(* --------------------------------------- *)
		EdMode == "Superc",
		QnsSectorList=Flatten[
			Outer[{##}&,##]&@@ConstantArray[
			Range[-L,L]
		,Norb],Max[1,Norb-1]],
(* --------------------------------------- *)
		EdMode == "InterorbSuperc",
		QnsSectorList = Range[-Norb*L,Norb*L];
	];
	QnsSectorList
];









End[];

EndPackage[];
