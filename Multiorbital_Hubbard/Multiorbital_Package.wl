(* ::Package:: *)

BeginPackage["DMFT`"]

StartingBath::usage = "StartingBath[InitializeBathMode, Nbath, Norb, EdMode] returns a list containing the bath parameters to start the DMFT loop.
If EdMode = ''Normal'' then the output has the form {e,V}, where e and V are lists of Norb x Nbath elements, representing the bath energies and the bath-impurity hybridizations.
If EdMode = ''Superc'' then the output has the form {e,V,\[CapitalDelta]}, where e and V are defined as above, and \[CapitalDelta] is the Norb x Nbath dimensional list of pairs creation (annihilation) amplitudes.
If EdMode = ''InterorbSuperc'' then the output has the form {e,V,\[CapitalDelta],\[CapitalXi]}, where e, V, \[CapitalDelta] are as above, and \[CapitalXi] is the Nbath - dimensional list of interorbital pairs creation (annihilation) amplitudes
InitializeBathMode is a string with the path to the file containing the bath parameters; if it is set to ''Default'', default parameters are dropped."

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


End[];

EndPackage[];
