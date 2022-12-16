(* ::Package:: *)

$Path = Join[$Path, {"C:\\Users\\matte\\Desktop\\Mathematica_package\\"}];
(* Import the code *)
<<("DMFT.wl");


(* Import the input file *)
<<("InputFile_Template.wl");
(* Prepare everything *)
<<("Preparation.wl");
BathParameters = {{{-1.203170050425068`,-0.12916886740504122`,-0.004458437824990661`,0.004457946089326225`,0.12916300033225653`,1.2031518006065318`},{-1.203170050425068`,-0.12916886740504122`,-0.004458437824990661`,0.004457946089326225`,0.12916300033225653`,1.2031518006065318`}},{{-0.29971575367547293`,0.14197278690924572`,0.042383222851114906`,0.04238173360852318`,0.14197073309327288`,-0.2997165371655205`},{-0.29971575367547293`,0.14197278690924572`,0.042383222851114906`,0.04238173360852318`,0.14197073309327288`,-0.2997165371655205`}}};
(* Define list of U values *)
Print["List of U values: "];
Ulist = Table[ConstantArray[U, Norb], {U, 4.2, 6.00, 0.2}]


Do[

	Print["------------------------------------------------------------------------------------"];
	Print[Style["Start the loop for U = " <> ToString[U], Blue, 24]];
	Print["------------------------------------------------------------------------------------"];

	(* Update interaction parameters *)
	If[HFMode, 
		shift = -U[[1]]/2 - (Ust + Usec)*(Norb - 1)/2.;,
	(* else *)
		shift = 0;
	];
	(* get flat list of interaction parameters *)
	InteractionParameters = Flatten[{\[Delta], U, Ust, Usec, Jph, Jse, - \[Mu] + shift}];
	(* Start DMFT Loops ... *)
	<<("DMFT_Loop.wl");

, {U, Ulist}]


ListPlot[{
	Im[\[CapitalSigma]]
	}, Joined->True, PlotRange->{{0,100},All}, PlotStyle->{Thick,Dashing[.05]}]


<<("Post_Processing.wl")


BathParameters


BathParameters
