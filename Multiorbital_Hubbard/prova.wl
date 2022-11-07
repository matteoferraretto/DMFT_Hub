(* ::Package:: *)

(* compute the Green function - <T [(c1* c_{up,orb1}(\[Tau]) + c2* cdg_{dw,orb2}(\[Tau])) (c1 cdg_{up,orb1}(0) + c2 c_{dw,orb2}(0)) ] > *)
GreenFunctionImpurity[L_, f_, Norb_, {orb1_,orb2_}, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_, OptionsPattern[] ] := Module[{
	(* compute cdg_{s,orb}|gs> and c_{s,orb}|gs> for s = up, dw *)
	adgup = ApplyCdg[L, f, Norb, 1, 1, orb1, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode],
    bdw = ApplyC[L, f, Norb, 1, 2, orb2, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode],
    aup = ApplyC[L, f, Norb, 1, 1, orb1, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode],
	bdgdw = ApplyCdg[L, f, Norb, 1, 2, orb2, Normalize[Gs], GsQns, Sectors, SectorsDispatch, EdMode],
    c1 = OptionValue[c1], c2 = OptionValue[c2],
	Odggs, Ogs, newqns, H, E0, a, b, GFOparticle, GFOhole, GFO},	
(* compute all the main contributions to the Green function *)
(*          G_O(z) "Particle" contribution             *)
	Odggs = c1*adgup + c2*bdw; (* apply Odg |gs> = (c1 * Cdg_{up,orb1} + c2 * C_{dw,orb2}) |gs> *)
	newqns = FinalSector[L, f, Norb, 1, 1, orb1, GsQns, "Creation", EdMode]; (* evaluate the quantum numbers of the final sector *)
	H = Hsectors[[newqns/.SectorsDispatch]]; (* Hamiltonian on that sector *)
	{E0,a,b} = Lanczos[H, Normalize[Odggs] ]; (* Apply Lanczos starting from Odg|gs> *)
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a] ]; (* Krylov matrix in the final sector *)
	GFOparticle = (Norm[Odggs]^2)*(
		InverseElement[
			SparseArray[(# + Egs) * IdentityMatrix[Length[a] ] - H]
		, {1, 1}] &/@ zlist);
(*           G_O(z) "Hole" contribution               *)
	Ogs = Conjugate[c1]*aup + Conjugate[c2]*bdgdw; (* apply O |gs> = (Conjugate[c1] * C_{up,orb1} + Conjugate[c2] * Cdg_{dw,orb2}) |gs> *)
	newqns = FinalSector[L, f, Norb, 1, 1, orb, GsQns, "Annihilation", EdMode]; (* evaluate the quantum numbers of the final sector *)
	H = Hsectors[[newqns/.SectorsDispatch]]; (* Hamiltonian on that sector *)
	{E0,a,b} = Lanczos[H, Normalize[Ogs] ]; (* Apply Lanczos starting from O|gs> *)
	H = SparseArray[DiagonalMatrix[b, 1] + DiagonalMatrix[b, -1] + DiagonalMatrix[a] ]; (* Krylov matrix in the final sector *)
	GFOhole = (Norm[Ogs]^2)*(
		InverseElement[
			SparseArray[(# - Egs) * IdentityMatrix[Length[a] ] + H]
		, {1, 1}] &/@ zlist);
	GFOparticle + GFOhole
];
Options[GreenFunctionImpurity] = {c1 -> 1., c2 -> 0.}; (* by default use Odg = adgup *)

GreenFunctionImpurityNambu[L_, f_, Norb_, Egs_, Gs_, GsQns_, Hsectors_, Sectors_, SectorsDispatch_, EdMode_, zlist_, OptionsPattern[] ] := Module[{
    NMatsubara = Length[zlist], zlistextended = Join[zlist, -zlist], orb = OptionValue[orb],
    GF, GFO, GFP
    },
    Which[ 
        EdMode == "Superc",
        (* initialize Green function as a NMatsubara x 2 x 2 tensor *)
        GF = ConstantArray[0, {NMatsubara, 2, 2}];
        (* if there is spin symmetry, then perform Lanczos only once, apply to +z and -z and then split the result *)
        {GF[[All,1,1]], GF[[All,2,2]]} = Partition[ 
            GreenFunctionImpurity[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
        NMatsubara]; 
        GF[[All,2,2]] = -GF[[All,2,2]]; (* GF_11 -> G_{up,orb1; up,orb1}(z) ;  GF_22 -> - G_{dw,orb1; dw,orb1}(-z) *) 
        (* compute GFO, where Odg = adg_up + a_dw *)
        GFO = GreenFunctionImpurity[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
        (* compute GFP, where Pdg = adg_up + I a_dw *)
        GFP = GreenFunctionImpurity[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
        (* compute the off-diagonal part using the diagonal part and GFO, GFP *)
	    GF[[All,1,2]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All,1,1]] + GF[[All,2,2]]));
	    GF[[All,2,1]] = (1./2.)*(GFO + I*GFP - 1.*(1 + I)*(GF[[All,1,1]] + GF[[All,2,2]]));
    (* ------------------------------------------------------------- *)
        EdMode == "InterorbSuperc",
        (* initialize Green function as a NMatsubara x 2 x 2 tensor *)
        GF = ConstantArray[0, {NMatsubara, 2*Norb, 2*Norb}];
        If[
        (* general case: NO orbital symmetry *)
            !OrbitalSymmetry,
            (*(* compute diagonal part of the diagonal blocks *)
            Do[
                {GF[[All, 2*(orb-1)+1, 2*(orb-1)+1]], GF[[All, 2*(orb-1)+2, 2*(orb-1)+2]]} = Partition[ 
                    GreenFunctionImpurity[L, f, Norb, {orb,orb}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
                NMatsubara]; 
                GF[[All, 2*(orb-1)+2, 2*(orb-1)+2]] = -GF[[All, 2*(orb-1)+2, 2*(orb-1)+2]]; 
            , {orb, Norb}];
            
            (* Speedup ? *)
            {GF[[All, 2*(#-1)+1, 2*(#-1)+1]], GF[[All, 2*(#-1)+2, 2*(#-1)+2]]} = Partition[ 
                    GreenFunctionImpurity[L, f, Norb, {#,#}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
                NMatsubara] &/@ Range[Norb];
            (GF[[All, 2*(#-1)+2, 2*(#-1)+2]] = -GF[[All, 2*(#-1)+2, 2*(#-1)+2]]) &/@ Range[Norb]; 
            *)
            (* compute the diagonal part of all the Nambu blocks *)
            Do[{
                GF[[All, 2(orb1-1)+1, 2(orb2-1)+1]], (* element 11, 13, 15, 33, 35, 55, ... *)
                GF[[All, 2(orb1-1)+2, 2(orb2-1)+2]] (* element 22, 24, 26, 44, 46, 66, ... *)
                } = Partition[ 
                    GreenFunctionImpurity[L, f, Norb, {orb1,orb2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
                NMatsubara]; 
                GF[[All, 2(orb1-1)+2, 2(orb2-1)+2]] = -GF[[All, 2(orb1-1)+2, 2(orb2-1)+2]]; (* put correct sign on elements 22, 24, 44, ...*)
            , {orb1, Norb}, {orb2, orb1, Norb}];
            Do[
				(* compute GFO, where Odg = cdg_{up,orb1} + c_{dw,orb2} *)
				GFO = GreenFunctionImpurity[L, f, Norb, {orb1,orb2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
				(* compute GFP, where Pdg = cdg_{up,orb1} + I c_{dw,orb2} *)
				GFP = GreenFunctionImpurity[L, f, Norb, {orb1,orb2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
				(* compute the top-right element of each block using the diagonal part and GFO, GFP *)
				GF[[All, 2(orb1-1)+1, 2(orb2-1)+2]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 2(orb1-1)+1, 2(orb1-1)+1]] + GF[[All, 2(orb2-1)+2, 2(orb2-1)+2]]));
				(* compute the bottom-left element of each non-diagonal block *)
				If[orb2 > orb1, 
					GFO = GreenFunctionImpurity[L, f, Norb, {orb2,orb1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
					GFP = GreenFunctionImpurity[L, f, Norb, {orb2,orb1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
					GF[[All, 2(orb1-1)+2, 2(orb2-1)+1]] = Conjugate[
                        (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 2(orb2-1)+1, 2(orb2-1)+1]] + GF[[All, 2(orb1-1)+2, 2(orb1-1)+2]]))
                    ]; (* you have to take the conjugate! For example GF23 is not Fba, but Fba* !  *)
				];
            , {orb1, Norb}, {orb2, orb1, Norb}];
            (* fill up the lower triangular part by conjugating the upper triangular part *)
            Do[
                GF[[All, i, j]] = Conjugate[GF[[All, j, i]]];
            , {i, 2Norb}, {j, i-1}];,
        (* ------------------------------------- *)
        (* else *)
            OrbitalSymmetry,
            (* compute the diagonal part of (a representative of) the diagonal Nambu block *)
            {GF[[All, 1, 1]], GF[[All, 2, 2]] } = 
            Partition[ 
                GreenFunctionImpurity[L, f, Norb, {1,1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
            NMatsubara]; 
            GF[[All, 2, 2]] = -GF[[All, 2, 2]];
            (* possibly GF13 = GF24 = 0 ? For now let's be safe ... *)
            (* compute the diagonal part of (a representative of) the diagonal Nambu block *)
            {GF[[All, 1, 3]], GF[[All, 2, 4]] } = 
            Partition[ 
                GreenFunctionImpurity[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlistextended],
            NMatsubara]; 
            GF[[All, 2, 4]] = -GF[[All, 2, 4]];
            (* *)
            (* compute GFO, where Odg = cdg_{up,1} + c_{dw,2} *)
			GFO = GreenFunctionImpurity[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
			(* compute GFP, where Pdg = cdg_{up,1} + I c_{dw,2} *)
			GFP = GreenFunctionImpurity[L, f, Norb, {1,2}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
            (* compute the top-right element of each block using the diagonal part and GFO, GFP *)
			GF[[All, 1, 2]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 1, 1]] + GF[[All, 2, 2]]));
            GF[[All, 1, 4]] = (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 1, 1]] + GF[[All, 2, 2]]));
            (* possibly GF23 = - GF14 ? For now let's be safe ... *)
            (* compute GFO, where Odg = cdg_{up,1} + c_{dw,2} *)
			GFO = GreenFunctionImpurity[L, f, Norb, {2,1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.];
			(* compute GFP, where Pdg = cdg_{up,1} + I c_{dw,2} *)
			GFP = GreenFunctionImpurity[L, f, Norb, {2,1}, Egs, Gs, GsQns, Hsectors, Sectors, SectorsDispatch, EdMode, zlist, c2 -> 1.*I];
            GF[[All, 2, 3]] = Conjugate[
                (1./2.)*(GFO - I*GFP - 1.*(1 - I)*(GF[[All, 1,1]] + GF[[All, 2, 2]]))
            ];
            (* fill up all the other entries according to the symmetries *)
            Do[
                GF[[All, 2(orb-1)+1, 2(orb-1)+1]] = GF[[All, 1, 1]];
                GF[[All, 2(orb-1)+2, 2(orb-1)+2]] = GF[[All, 2, 2]];
            , {orb, 2, Norb}]; (* diagonal elements of diagonal blocks *)
            Do[
                If[orb1 == 1 && orb2 == 2, Continue[]; (* this is done already! *) ];
                GF[[All, 2(orb1-1)+1, 2(orb2-1)+1]] = GF[[All, 1, 3]];
                GF[[All, 2(orb1-1)+2, 2(orb2-1)+2]] = GF[[All, 2, 4]];
                GF[[All, 2(orb1-1)+1, 2(orb2-1)+2]] = GF[[All, 1, 4]];
                GF[[All, 2(orb1-1)+2, 2(orb2-1)+3]] = GF[[All, 2, 3]];
            , {orb1, Norb}, {orb2, orb1+1, Norb}]; (* elements of non-diagonal blocks *)
            Do[
                GF[[All, i, j]] = Conjugate[GF[[All, j, i]]];
            , {i, 2Norb}, {j, i-1}]; (* fill up the lower triangular part by conjugating the upper triangular part *)
        ]
    ];
    GF
];
Options[GreenFunctionImpurityNambu] = {orb -> 1}; (* if orbitals are not symmetric and EdMode = "Superc", you can specify the orbital index *)

(*
    to run this script:
        - from the terminal navigate to the file directory
        - type
            $ wolframscript -file prova.wl
        - in the code, explicitly use Print[] to show output
*)
