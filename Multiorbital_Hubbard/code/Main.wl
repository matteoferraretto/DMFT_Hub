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


<<("Post_Processing.wl");


ListPlot[spectralfunction[[1]],PlotRange->{{-6,4},All},Joined->True,AspectRatio->1/3]


2*Sum[spectralfunction[[2]][[All,2]][[j]], {j, 5005}]*d\[Omega]
2*Sum[spectralfunction[[2]][[All,2]][[j]], {j, 5000, 10000}]*d\[Omega]
2*(1.-0.9986922573937804`)


ListPlot[
{\[Omega],Im[\[CapitalSigma]realfreq[[1]]]}\[Transpose], Joined->True, PlotRange->All
]


spec = -(1./Pi) * Im[Mean[MapApply[
	GreenFunctionImpurity[L, f, Norb, {1,1}, Egs, ##, Hsectors, Sectors, SectorsDispatch, EdMode, \[Omega]+I*\[Eta]]&,
	{Gs, GsQns}\[Transpose]
]]];

spec0 = - (1./Pi) * Im[1./((
	(Weiss/.{\[Mu]eff -> \[Mu] - \[Delta][[1]]})/.Thread[symbols -> TakeIndependentParameters[L, f, Norb, 1, 1, BathParameters, EdMode]])/.{z -> #}&/@(\[Omega] + I*\[Eta]))
];

ListPlot[{\[Omega],spec}\[Transpose], Joined->True, PlotRange->All]
ListPlot[{\[Omega],spec0}\[Transpose], Joined->True, PlotRange->All]

Total[spec]*d\[Omega]
Total[spec0]*d\[Omega]


latticespec = -(1./Pi) * Sum[ 
	LatticeWeights[[k]]*( - \[Eta])/((\[Omega] + \[Mu] - LatticeEnergies[[k,1,1]] )^2 + \[Eta]^2)
, {k, LatticePoints}];

ListPlot[{\[Omega],latticespec}\[Transpose], Joined->True, PlotRange->All]


Total[latticespec]*d\[Omega]


ListPlot[Im[\[CapitalSigma][[1]]], PlotRange->All, Joined->True]


ListPlot[Im[\[CapitalSigma]], Joined->True, PlotRange->{{0,500},All}]


\[CapitalSigma]Adriano = Import["C:\\Users\\matte\\Desktop\\Ph.D\\PhD Project\\Assignment1\\Hubbard-model\\U3.40\\impSigma_l11_s1_realw.ed", "Table"];

loc = Sum[
	LatticeWeights[[i]]/(\[CapitalSigma]Adriano[[All,1]] + \[Mu] - LatticeEnergies[[All,1,1]][[i]] - \[CapitalSigma]Adriano[[All,3]] - I*\[CapitalSigma]Adriano[[All,2]])
, {i, LatticePoints}];

ListPlot[{
	\[CapitalSigma]Adriano[[All,1]],
	-(1./Pi)*Im[loc]
	}\[Transpose]
, Joined->True, PlotRange->All]



grid =  Tuples[{
		Range[ 1500 ],
		Range[ 10000 ]
	}];
	
Dimensions[Partition[grid, 1500]]
