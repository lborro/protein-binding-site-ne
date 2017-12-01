package br.embrapa.cnptia.gpbc.plc.descriptors.plp;

import org.biojava.nbio.structure.Element;

import br.embrapa.cnptia.cbi.sdl.core.Bond;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;

public class AtomTypingRules {
	public static final int ACCEPTOR = 1;
	public static final int DONOR = 2;
	public static final int ACCEPTOR_DONOR = 3;
	public static final int METAL = 4;
	public static final int NON_POLAR = 5;

	public static int getAtomType(IAtom atom) {

		if(atom.getElement().equals(Element.N) || atom.getElement().equals(Element.O)){
			int nHydrogens = 0;
			for (Bond bond : atom.getBonds()){
				IAtom other = bond.getOther(atom);
				if(other.getElement().equals(Element.H)) nHydrogens++; 
			}
			if(nHydrogens == 0) return ACCEPTOR;
			if(atom.getElement().equals(Element.N)) return DONOR;
			if(nHydrogens == 1 || atom.getResidue().isWater()) return ACCEPTOR_DONOR;
		}
		if(atom.getElement().isMetal() ) return METAL;

		return NON_POLAR;
	}
}
