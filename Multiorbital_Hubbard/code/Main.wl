(* ::Package:: *)

$Path = Join[$Path, {"C:\\Users\\matte\\Desktop\\Mathematica_package\\"}];
(* Import the code *)
<<"DMFT.wl";


(* Import the input file *)
<<"InputFile_Template.wl";


(* Prepare everything *)
<<"Preparation.wl";


(* Start DMFT Loops ... *)
Off[Eigensystem::arh];
<<"DMFT_Loop.wl";


<<"Post_Processing.wl";
