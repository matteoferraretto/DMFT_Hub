(* ::Package:: *)

$Path = Join[$Path, {"C:\\Users\\matte\\Desktop\\Mathematica_package\\"}];
(* Import the code *)
<<"DMFT.wl";


(* Import the input file *)
<<"InputFile_Template.wl";


(* Prepare everything *)
<<"Preparation.wl";


(* Start DMFT Loops ... *)
<<"DMFT_Loop.wl";


(*<<"Post_Processing.wl";*)
<<"Observables.wl";
(*KineticEnergyNew[\[Mu], LatticeEnergies, LatticeWeights, \[CapitalSigma], i\[Omega], EdMode]*)


(* \[CapitalSigma]0 = Last /@ \[CapitalSigma];

2. * 2. * TMats * Total @ Table[
	LatticeWeights[[k]] * Tr[ (LatticeEnergies[[k]] - \[Mu]) * Total @ Table[
		Inverse[(i\[Omega][[n]] + \[Mu])*IdentityMatrix[Norb] - LatticeEnergies[[k]] - DiagonalMatrix[\[CapitalSigma][[All, n]]]]
		- (LatticeEnergies[[k]] + \[CapitalSigma]0)/(i\[Omega][[n]]^2)
		, {n, CGNMatsubara}]]
	, {k, LatticePoints}] + 
	0.5 * Sum[
		LatticeWeights[[k]] * Tr[(LatticeEnergies[[k]] - \[Mu])]
	, {k, LatticePoints}] - 
	(1./(2.*TMats)) * Sum[
		LatticeWeights[[k]] * Tr[(LatticeEnergies[[k]] - \[Mu]) * (LatticeEnergies[[k]] + \[CapitalSigma]0)]
	, {k, LatticePoints}] *)


(* non interacting kinetic energy *)
N[ 2*2*Integrate[(2/(Pi))*x*Sqrt[1-x^2], {x, -1, 0}] ]


ListPlot[{
	Re[ \[CapitalSigma][[1,All,1,1]] ],
	Re[ \[CapitalSigma][[1,All,2,2]] ]
}, Joined->True, PlotStyle->{Thick,Dashing[.05]}]



HlocBlocks[[13]]//Dimensions


Length[Sectors]
Length[Flatten[InteractionParameters]]



NewBathParameters = BathParameters;


(* {\[Delta], \[Mu](n=2.0)} *)
(* {{0, 0}, {0.15, 0.075}} *)


udg = (Normalize[#] &/@
Eigenvectors[
PauliMatrix[1]
]\[Transpose]);
u = ConjugateTranspose[udg];

udg // MatrixForm
u // MatrixForm

u . PauliMatrix[1] . udg

FullSimplify[
(u . {
{a, 0},
{0, b}
} . udg)] // MatrixForm
