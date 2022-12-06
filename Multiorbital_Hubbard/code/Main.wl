(* ::Package:: *)

$Path = Join[$Path, {"C:\\Users\\matte\\Desktop\\Mathematica_package\\"}];
(* Import the code *)
<<("DMFT.wl");


(* Import the input file *)
<<("InputFile_Template.wl");


(* Prepare everything *)
<<("Preparation.wl");
BathParameters = {{{-0.014254812983406195`,-0.000486000389943534`,0.000264730758641006`,0.012798715205599743`},{-0.014254812983406195`,-0.000486000389943534`,0.000264730758641006`,0.012798715205599743`}},{{0.31502007578654906`,0.11280704199916472`,0.10485550277818963`,0.3176257055691743`},{0.31502007578654906`,0.11280704199916472`,0.10485550277818963`,0.3176257055691743`}},{{-0.33497341228945526`,-0.026217082713124513`,0.017885553839128123`,0.30891934932789067`}}};


(* Start DMFT Loops ... *)
<<("DMFT_Loop.wl");


<<("Post_Processing.wl");


ListDensityPlot[
	- (1./Pi) * Im[
		TwoByTwoInverse[
			Table[
				(\[Omega][[i]]+I*\[Eta])*IdentityMatrix[2] + (\[Mu] - LatticeEnergies[[All,1,1]][[Range[1500]]][[j]])*PauliMatrix[3] - \[CapitalSigma]realfreq[[1]][[i]]
			,{i, 1, 10000,20}, {j, 1, 1500}]
		][[All,All,1,1]]
	]
]
