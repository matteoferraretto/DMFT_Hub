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
(*
messageHandler = If[Last[#], Abort[]] &;
Internal`AddHandler["Message", messageHandler]
*)
<<"DMFT_Loop.wl";


(*<<"Post_Processing.wl";*)


(1.6970055054809476`- 0.3029949986056395`)/2.


0.036972233810872654` -1.9630266665044207`
