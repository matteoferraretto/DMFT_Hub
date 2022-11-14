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


L = 3; f = 2; Norb = 2; EdMode = "FullSuperc"; \[Mu] = 0;

symbols = Symbols[L, f, Norb, EdMode]

Apart[WeissField[L, f, Norb, \[Mu], symbols, z, EdMode],z]
