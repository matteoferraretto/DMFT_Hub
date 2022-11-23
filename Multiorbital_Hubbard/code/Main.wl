(* ::Package:: *)

$Path = Join[$Path, {"C:\\Users\\matte\\Desktop\\Mathematica_package\\"}];
(* Import the code *)
<<("DMFT.wl");


(* Import the input file *)
<<("InputFile_Template.wl");


(* Prepare everything *)
<<("Preparation.wl");


(* Start DMFT Loops ... *)
<<("DMFT_Loop.wl");


<<("Post_Processing.wl")


Flatten @ Table[
	Print[\[Alpha]," ",\[Beta]];
	Total[
		Abs[InverseG0[[All, \[Alpha], \[Beta]]] - InverseG0old[[All, \[Alpha], \[Beta]]]]
	]/Max[
		Total[Abs[InverseG0[[All,\[Alpha],\[Beta]]]]], Total[Abs[InverseG0old[[All,\[Alpha],\[Beta]]]]]
	]
, {\[Alpha], 4}, {\[Beta], \[Alpha], 4}]

Mean[%]

Total[Abs[InverseG0[[All,1,3]]]]
Total[Abs[InverseG0old[[All,1,3]]]]
Total[Abs[InverseG0[[All,1,3]]] - Abs[InverseG0old[[All,1,3]]]]


Weiss/.{\[CapitalXi]1->0, \[CapitalXi]2->0}


MapApply[
	InverseGreenFunction[L, f, Norb, 1, 1, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, i\[Omega]]&,
{Gs, GsQns}\[Transpose]]

