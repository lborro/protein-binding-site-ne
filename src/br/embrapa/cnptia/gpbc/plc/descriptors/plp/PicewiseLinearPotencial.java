package br.embrapa.cnptia.gpbc.plc.descriptors.plp;

import java.util.List;

import br.embrapa.cnptia.gpbc.plc.structure.ProteinLigandAtomPair;

public class PicewiseLinearPotencial {
	public static final double HB_A =  2.3;
	public static final double HB_B =  2.6;
	public static final double HB_C =  3.1;
	public static final double HB_D =  3.4;
	public static final double HB_E = -1.0;
	public static final double HB_F = 20.0;

	public static final double METAL_A =  1.4;
	public static final double METAL_B =  2.2;
	public static final double METAL_C =  2.6;
	public static final double METAL_D =  2.8;
	public static final double METAL_E = -1.0;
	public static final double METAL_F = 20.0;

	public static final double BURIED_A =  3.4;
	public static final double BURIED_B =  3.6;
	public static final double BURIED_C =  4.5;
	public static final double BURIED_D =  5.5;
	public static final double BURIED_E = -0.1;
	public static final double BURIED_F = 20.0;

	public static final double NON_POLAR_A =  3.4;
	public static final double NON_POLAR_B =  3.6;
	public static final double NON_POLAR_C =  4.5;
	public static final double NON_POLAR_D =  5.5;
	public static final double NON_POLAR_E = -0.4;
	public static final double NON_POLAR_F = 20.0;

	public static final double REP_A =  3.2;
	public static final double REP_B =  5.0;
	public static final double REP_C =  0.1;
	public static final double REP_D = 20.0;


	private static double getPLP(double r, double A, double B, double C, double D, double E, double F){
		if(r < A)
			return F*(A - r)/A;
		else if(r < B)
			return E*(r-A)/(B-A);
		else if(r < C)
			return E;
		else if(r <= D) 
			return E*(D-r)/(D-C);		
		return 0;
	}

	private static double getRepulsion(double r, double A, double B, double C, double D) {
		if(r < A)
			return r*(C-D)/A + D;
		else if (r < B) 
			return -C*(r-A)/(B-A) + C;
		return 0;
	}

	public static double calculate(List<ProteinLigandAtomPair> pairs) {
		double plp = 0;
		for(ProteinLigandAtomPair pair: pairs){
			int interactionType = InteractionRules.getAtomAtomPairInteraction(pair);

			//Hydrogen bond 
			if(interactionType == InteractionRules.HYDROGEN_BOND){
				plp += PicewiseLinearPotencial.getPLP(
						pair.getDistance(),
						PicewiseLinearPotencial.HB_A,
						PicewiseLinearPotencial.HB_B,
						PicewiseLinearPotencial.HB_C,
						PicewiseLinearPotencial.HB_D,
						PicewiseLinearPotencial.HB_E,
						PicewiseLinearPotencial.HB_F);
				//Metal interaction
			} else if(interactionType == InteractionRules.METAL){
				plp += PicewiseLinearPotencial.getPLP(
						pair.getDistance(),
						PicewiseLinearPotencial.METAL_A,
						PicewiseLinearPotencial.METAL_B,
						PicewiseLinearPotencial.METAL_C,
						PicewiseLinearPotencial.METAL_D,
						PicewiseLinearPotencial.METAL_E,
						PicewiseLinearPotencial.METAL_F);

				//Non-polar interaction
			} else if(interactionType == InteractionRules.NON_POLAR){
				plp += PicewiseLinearPotencial.getPLP(
						pair.getDistance(),
						PicewiseLinearPotencial.NON_POLAR_A,
						PicewiseLinearPotencial.NON_POLAR_B,
						PicewiseLinearPotencial.NON_POLAR_C,
						PicewiseLinearPotencial.NON_POLAR_D,
						PicewiseLinearPotencial.NON_POLAR_E,
						PicewiseLinearPotencial.NON_POLAR_F);
				//Buried interaction
			} else if(interactionType == InteractionRules.BURIED) {
				plp += PicewiseLinearPotencial.getPLP(
						pair.getDistance(),
						PicewiseLinearPotencial.BURIED_A,
						PicewiseLinearPotencial.BURIED_B,
						PicewiseLinearPotencial.BURIED_C,
						PicewiseLinearPotencial.BURIED_D,
						PicewiseLinearPotencial.BURIED_E,
						PicewiseLinearPotencial.BURIED_F);
				//Repulsion
			} else {
				plp += PicewiseLinearPotencial.getRepulsion(
						pair.getDistance(),
						PicewiseLinearPotencial.REP_A,
						PicewiseLinearPotencial.REP_B,
						PicewiseLinearPotencial.REP_C,
						PicewiseLinearPotencial.REP_D);
			}
		}
		return plp;
	}
}
