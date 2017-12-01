package br.embrapa.cnptia.gpbc.plc.descriptors.plp;

import br.embrapa.cnptia.gpbc.plc.structure.ProteinLigandAtomPair;

public class InteractionRules {
	public static final int HYDROGEN_BOND = 1;
	public static final int REPULSIVE = 2;
	public static final int BURIED = 3;
	public static final int METAL = 4;
	public static final int NON_POLAR = 5; 

	public static int getAtomAtomPairInteraction(ProteinLigandAtomPair pair) {
		int ligAtomType = AtomTypingRules.getAtomType(pair.getLigandAtom());
		int proAtomType = AtomTypingRules.getAtomType(pair.getProteinAtom());

		int type;

		//Ligand Atom: Donor
		if(ligAtomType == AtomTypingRules.DONOR) {
			if(proAtomType == AtomTypingRules.DONOR || proAtomType ==  AtomTypingRules.METAL)
				type = InteractionRules.REPULSIVE;
			else if(proAtomType == AtomTypingRules.ACCEPTOR || proAtomType == AtomTypingRules.ACCEPTOR_DONOR)
				type = InteractionRules.HYDROGEN_BOND;
			else 
				type = InteractionRules.BURIED;
		} 
		//Ligand Atom: Acceptor
		else if (ligAtomType == AtomTypingRules.ACCEPTOR){
			if(proAtomType == AtomTypingRules.DONOR || proAtomType == AtomTypingRules.ACCEPTOR_DONOR)
				type = InteractionRules.HYDROGEN_BOND;
			else if(proAtomType == AtomTypingRules.ACCEPTOR)
				type = InteractionRules.REPULSIVE;
			else if(proAtomType == AtomTypingRules.METAL)
				type = InteractionRules.METAL;
			else
				type = InteractionRules.BURIED;
		} 
		//Ligand Atom: Acceptor & Donor
		else if(ligAtomType == AtomTypingRules.ACCEPTOR_DONOR) {
			if(proAtomType == AtomTypingRules.METAL)
				type = InteractionRules.METAL;
			else if(proAtomType == AtomTypingRules.NON_POLAR)
				type = InteractionRules.BURIED;
			else 
				type = InteractionRules.HYDROGEN_BOND;
		} 
		//Ligand Atom: Non polar
		else {
			if(proAtomType == AtomTypingRules.NON_POLAR)
				type = InteractionRules.NON_POLAR;
			else
				type = InteractionRules.BURIED;
		}
		return type;
	}
}
